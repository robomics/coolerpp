// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_test_macros.hpp>

#include "coolerpp/benchmark/common.hpp"
#include "coolerpp/chromosome.hpp"

namespace coolerpp::benchmark {

TEST_CASE("ChromosomeSet: benchmark", "[chroms][bench]") {
  const auto chroms = constants::hg38_chroms;

  BENCHMARK_ADVANCED("ChromosomeSet: contains(Chromosome)")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
    const auto chrom_list = generate_random_chrom_list(chroms, num_runs);
    meter.measure([&](std::size_t i) { return chroms.contains(chrom_list[i]); });
  };

  BENCHMARK_ADVANCED("ChromosomeSet: at(id)")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
    const auto chrom_ids = generate_random_chrom_id_list(chroms, num_runs);
    meter.measure([&](std::size_t i) { return chroms.at(chrom_ids[i]); });
  };

  BENCHMARK_ADVANCED("ChromosomeSet: at(name)")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
    const auto chrom_list = generate_random_chrom_list(chroms, num_runs);
    meter.measure([&](std::size_t i) { return chroms.at(chrom_list[i].name()); });
  };

  BENCHMARK_ADVANCED("ChromosomeSet: get_id(Chromosome)")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
    const auto chrom_list = generate_random_chrom_list(chroms, num_runs);
    meter.measure([&](std::size_t i) { return chroms.get_id(chrom_list[i].name()); });
  };
}

}  // namespace coolerpp::benchmark
