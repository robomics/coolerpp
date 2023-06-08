// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <random>
#include <string_view>
#include <vector>

#include "../../../../../test/units/include/coolerpp/test/self_deleting_folder.hpp"
#include "coolerpp/chromosome.hpp"

namespace coolerpp::benchmark {

namespace constants {
inline constexpr std::array<std::string_view, 24> chrom_names{
    {"chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",
     "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
     "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",  "chrY"}};

inline constexpr std::array<std::uint32_t, 24> chrom_sizes{
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
    138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
    83257441,  80373285,  58617616,  64444167,  46709983,  50818468,  156040895, 57227415,
};

// NOLINTNEXTLINE
inline const ChromosomeSet hg38_chroms{chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()};

}  // namespace constants

using SelfDeletingFolder = test::SelfDeletingFolder;

template <typename It>
[[nodiscard]] inline std::mt19937_64 get_prng(It first_seed, It last_seed) {
  auto seq = std::seed_seq(first_seed, last_seed);
  return std::mt19937_64{seq};
}

template <typename I>
[[nodiscard]] inline std::mt19937_64 get_prng(const std::vector<I>& seeds) {
  return get_prng(seeds.begin(), seeds.end());
}

[[nodiscard]] inline std::mt19937_64 get_prng() {
  std::random_device rd{};
  auto seq = std::seed_seq({rd(), rd(), rd(), rd(), rd(), rd()});
  return std::mt19937_64{seq};
}

template <typename I>
[[nodiscard]] inline std::mt19937_64 get_prng(I seed) {
  static_assert(std::is_integral_v<I>);
  auto seq = std::seed_seq(seed);
  return std::mt19937_64{seq};
}

[[nodiscard]] inline std::vector<std::uint32_t> generate_random_chrom_id_list(
    const ChromosomeSet& chroms, std::size_t size) {
  std::vector<std::uint32_t> chrom_sizes(chroms.size());
  std::transform(chroms.begin(), chroms.end(), chrom_sizes.begin(),
                 [](const Chromosome& c) { return c.size(); });

  auto rand_eng = get_prng(chrom_sizes.begin(), chrom_sizes.end());

  std::vector<std::uint32_t> chrom_ids(size);
  std::generate(chrom_ids.begin(), chrom_ids.end(), [&]() {
    const std::uint32_t lb = 0;
    const auto ub = static_cast<std::uint32_t>(chroms.size() - 1);
    return std::uniform_int_distribution<std::uint32_t>{lb, ub}(rand_eng);
  });

  return chrom_ids;
}

[[nodiscard]] inline std::vector<Chromosome> generate_random_chrom_list(const ChromosomeSet& chroms,
                                                                        std::size_t size) {
  const auto chrom_ids = generate_random_chrom_id_list(chroms, size);

  std::vector<Chromosome> chrom_list(size);
  std::transform(chrom_ids.begin(), chrom_ids.end(), chrom_list.begin(),
                 [&](auto id) { return chroms.at(id); });

  return chrom_list;
}

}  // namespace coolerpp::benchmark
