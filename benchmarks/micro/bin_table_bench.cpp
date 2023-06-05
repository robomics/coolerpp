// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <array>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_test_macros.hpp>

#include "coolerpp/benchmark/common.hpp"
#include "coolerpp/bin_table.hpp"

namespace coolerpp::benchmark {

TEST_CASE("BinTable: benchmark", "[bins][bench]") {
  const auto chroms = constants::hg38_chroms;
  constexpr std::array<std::uint32_t, 8> resolutions{1,      10,      100,       1'000,
                                                     10'000, 100'000, 1'000'000, 10'000'000};

  for (const auto res : resolutions) {
    BENCHMARK_ADVANCED("BinTable: ctor(ChromosomeSet, bin_size) " + std::to_string(res))
    (Catch::Benchmark::Chronometer meter) {
      const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
      std::vector<Catch::Benchmark::storage_for<BinTable>> storage(num_runs);
      meter.measure([&](std::size_t i) { storage[i].construct(chroms, res); });
    };

    BENCHMARK_ADVANCED("BinTable: subset(Chromosome) " + std::to_string(res))
    (Catch::Benchmark::Chronometer meter) {
      const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
      const BinTable bins(chroms, res);

      const auto chrom_list = generate_random_chrom_list(chroms, num_runs);
      meter.measure([&](std::size_t i) { return bins.subset(chrom_list[i]); });
    };

    BENCHMARK_ADVANCED("BinTable: at(id) " + std::to_string(res))
    (Catch::Benchmark::Chronometer meter) {
      const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
      const BinTable bins(chroms, res);

      const auto chrom_ids = generate_random_chrom_id_list(chroms, num_runs);
      meter.measure([&](std::size_t i) { return bins.at(chrom_ids[i]); });
    };
    BENCHMARK_ADVANCED("BinTable: at_hint(id, Chromosome) " + std::to_string(res))
    (Catch::Benchmark::Chronometer meter) {
      const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
      const BinTable bins(chroms, res);

      const auto chrom_list = generate_random_chrom_list(bins.chromosomes(), num_runs);
      std::vector<std::uint64_t> bin_ids(num_runs);

      std::transform(chrom_list.begin(), chrom_list.end(), bin_ids.begin(),
                     [&, rand_eng = get_prng()](const Chromosome& chrom) mutable {
                       const auto pos = std::uniform_int_distribution<std::uint32_t>{
                           0, chrom.size() - 1}(rand_eng);
                       return bins.at(chrom, pos).id();
                     });

      meter.measure([&](std::size_t i) { return bins.at_hint(bin_ids[i], chrom_list[i]); });
    };
  }
}

}  // namespace coolerpp::benchmark
