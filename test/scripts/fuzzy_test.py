#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import ctypes
import itertools
import logging
import multiprocessing as mp
import pathlib
import random
import shlex
import shutil
import subprocess as sp
import sys
import time
from typing import Dict, Tuple

import cooler
import numpy as np
import pandas as pd


def make_cli():
    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    def valid_fraction(arg):
        if (n := float(arg)) >= 0 and n <= 1:
            return n

        raise ValueError("Not a number between 0 and 1")

    def valid_executable(arg):
        if (cmd := shutil.which(arg)) is not None:
            return pathlib.Path(cmd)

        raise FileNotFoundError(f'Unable to find executable "{arg}"')

    cli = argparse.ArgumentParser()

    cli.add_argument(
        "cooler",
        type=pathlib.Path,
        help="Path to a Cooler file (URI syntax supported).",
    )

    cli.add_argument(
        "--1d-to-2d-query-ratio",
        type=valid_fraction,
        default=0.33,
        help="Ratio of 1D to 2D queries. Use 0 or 1 to only test 1D or 2D queries.",
    )

    cli.add_argument("--duration", type=positive_int, default=60, help="Duration in seconds.")
    cli.add_argument(
        "--path-to-coolerpp-dump", type=valid_executable, default="coolerpp_dump", help="Path to coolerpp_dump binary."
    )
    cli.add_argument("--query-length-avg", type=float, default=5_000_000, help="Average query size.")
    cli.add_argument("--query-length-std", type=float, default=1_000_000, help="Standard deviation for query size.")
    cli.add_argument("--seed", type=int, default=2074288341)
    cli.add_argument(
        "--nproc",
        type=int,
        choices=range(1, mp.cpu_count() + 1),
        default=mp.cpu_count(),
        help="Maximum number of parallel processes.",
    )
    return cli


def cooler_dump(selector, query1: str, query2: str):
    logging.debug("[cooler] running query for %s, %s...", query1, query2)
    return selector.fetch(query1, query2)


def coolerpp_dump(coolerpp_bin: pathlib.Path, path_to_cooler_file: pathlib.Path, query1: str, query2: str):
    cmd = [shutil.which(str(coolerpp_bin)), str(path_to_cooler_file), query1, query2]

    cmd = shlex.split(" ".join(str(tok) for tok in cmd))
    logging.debug("[coolerpp] Running %s...", cmd)
    sp.check_output(cmd, stderr=sp.STDOUT)

    with sp.Popen(cmd, stdin=None, stderr=sp.PIPE, stdout=sp.PIPE) as coolerpp_dump:
        df = pd.read_table(
            coolerpp_dump.stdout, names=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"]
        )
        coolerpp_dump.communicate()
        if (code := coolerpp_dump.returncode) != 0:
            print(coolerpp_dump.stderr, file=sys.stderr)
            raise RuntimeError(f"{cmd} terminated with code {code}")

    return df


def read_chrom_sizes(path_to_cooler_file: pathlib.Path) -> Dict[str, int]:
    return cooler.Cooler(str(path_to_cooler_file)).chromsizes.to_dict()


def generate_query_1d(chroms, weights: np.ndarray, mean_length: float, stddev_length: float) -> str:
    chrom_name, chrom_size = random.choices(chroms, weights=weights, k=1)[0]

    query_length = max(2.0, random.gauss(mu=mean_length, sigma=stddev_length))

    center_pos = random.randint(0, chrom_size)
    start_pos = max(0.0, center_pos - (query_length / 2))
    end_pos = min(chrom_size, start_pos + query_length)

    return f"{chrom_name}:{start_pos:.0f}-{end_pos:.0f}"


def generate_query_2d(
    chroms, weights: np.ndarray, ranks: Dict[str, int], mean_length: float, stddev_length: float
) -> Tuple[str, str]:
    q1 = generate_query_1d(chroms, weights, mean_length, stddev_length)
    q2 = generate_query_1d(chroms, weights, mean_length, stddev_length)

    chrom1, _, coord1 = q1.partition(":")
    chrom2, _, coord2 = q2.partition(":")

    if ranks[chrom1] > ranks[chrom2]:
        q1, q2 = q2, q1

    if chrom1 == chrom2:
        start1, _, _ = coord1.partition("-")
        start2, _, _ = coord2.partition("-")
        if int(start1) > int(start2):
            q1, q2 = q2, q1

    return q1, q2


def results_are_identical(worker_id, q1, q2, expected, found) -> bool:
    if len(expected) != len(found):
        logging.warning(
            "[%d] %s, %s: FAIL! Expected %d nnz, found %d!",
            worker_id,
            q1,
            q2,
            len(expected),
            len(found),
        )
        return False

    diff = expected.compare(found)

    if len(diff) != 0:
        logging.warning(
            "[%d] %s, %s (%d nnz): FAIL! Found %d differences!",
            worker_id,
            q1,
            q2,
            len(expected),
            len(diff),
        )
        return False

    logging.debug("[%d] %s, %s (%d nnz): OK!", worker_id, q1, q2, len(expected))
    return True


def seed_prng(worker_id: int, seed):
    seed = hash(tuple([worker_id, seed]))
    logging.info("[%d] seed: %d", worker_id, seed)
    random.seed(seed)


def worker(
    path_to_cooler: pathlib.Path,
    path_to_coolerpp_dump: pathlib.Path,
    chroms_flat,
    chrom_ranks,
    query_length_mu: float,
    query_length_std: float,
    _1d_to_2d_query_ratio: float,
    seed: int,
    worker_id: int,
    end_time,
) -> Tuple[int, int]:
    global early_return

    num_failures = 0
    num_queries = 0

    try:
        seed_prng(worker_id, seed)

        chrom_sizes = np.array([n for _, n in chroms_flat], dtype=int)
        weights = chrom_sizes / chrom_sizes.sum()

        sel = cooler.Cooler(str(path_to_cooler)).matrix(balance=False, as_pixels=True, join=True)

        while time.time() < end_time:
            if early_return.value:
                logging.debug("[%d] early return signal received. Returning immediately!", worker_id)
                break

            if _1d_to_2d_query_ratio <= random.random():
                q1, q2 = generate_query_2d(
                    chroms_flat, weights, chrom_ranks, mean_length=query_length_mu, stddev_length=query_length_std
                )
            else:
                q1 = generate_query_1d(
                    chroms_flat, weights, mean_length=query_length_mu, stddev_length=query_length_std
                )
                q2 = q1

            num_queries += 1
            expected = cooler_dump(sel, q1, q2)
            found = coolerpp_dump(path_to_coolerpp_dump, path_to_cooler, q1, q2)

            if not results_are_identical(worker_id, q1, q2, expected, found):
                num_failures += 1

    except:
        logging.debug("[%d] exception raised in worker process. Sending early return signal!", worker_id)
        early_return.value = True
        raise

    return num_queries, num_failures


def main():
    args = vars(make_cli().parse_args())

    chroms = read_chrom_sizes(args["cooler"])
    chrom_ranks = {chrom: i for i, chrom in enumerate(chroms.keys())}
    chroms_flat = list(chroms.items())

    end_time = time.time() + args["duration"]

    with mp.Pool(args["nproc"]) as pool:
        results = pool.starmap(
            worker,
            zip(
                itertools.repeat(args["cooler"]),
                itertools.repeat(args["path_to_coolerpp_dump"]),
                itertools.repeat(chroms_flat),
                itertools.repeat(chrom_ranks),
                itertools.repeat(args["query_length_avg"]),
                itertools.repeat(args["query_length_std"]),
                itertools.repeat(args["1d_to_2d_query_ratio"]),
                itertools.repeat(args["seed"]),
                range(args["nproc"]),
                itertools.repeat(end_time),
            ),
            chunksize=1,
        )

    num_queries = sum((n for n, _ in results))
    num_failures = sum((n for _, n in results))
    num_passes = num_queries - num_failures
    if num_failures == 0:
        lvl = logging.INFO
    else:
        lvl = logging.WARN

    logging.log(
        lvl, "Score: %.4g%% (%d successes and %d failures).", 100 * num_passes / num_queries, num_passes, num_failures
    )

    return num_failures != 0


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    setup_logger()
    early_return = mp.Value(ctypes.c_bool, False)
    sys.exit(main())
