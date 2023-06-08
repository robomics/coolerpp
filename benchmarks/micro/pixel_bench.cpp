// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <array>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "coolerpp/benchmark/common.hpp"
#include "coolerpp/coolerpp.hpp"
#include "coolerpp/pixel.hpp"

namespace coolerpp::benchmark {

inline const std::filesystem::path datadir{"benchmarks/data"};  // NOLINT(cert-err58-cpp)

template <typename N>
[[nodiscard]] static std::vector<Pixel<N>> random_sample_pixels(
    std::mt19937_64& rand_eng, const std::vector<Pixel<N>>& population, std::size_t sample_size) {
  std::vector<Pixel<N>> buffer(sample_size);

  // We can't use std::sample because sample_size could be > population
  std::generate(buffer.begin(), buffer.end(), [&]() {
    const auto i = std::uniform_int_distribution<std::size_t>{0, population.size() - 1}(rand_eng);
    return population[i];
  });
  return buffer;
}

TEST_CASE("Pixel<N>: benchmark", "[pixel][bench]") {
  const auto test_file = datadir / "4DNFI9FVHJZQ.0.9.1.mcool::/resolutions/1000000";

  const auto clr = File::open_read_only(test_file.string());
  const std::vector<Pixel<std::int32_t>> pixels_pop{clr.begin<std::int32_t>(),
                                                    clr.end<std::int32_t>()};

  BENCHMARK_ADVANCED("Pixel<N>: operator==")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
    auto rand_eng = get_prng();
    const auto pixels1 = random_sample_pixels(rand_eng, pixels_pop, num_runs);
    const auto pixels2 = random_sample_pixels(rand_eng, pixels_pop, num_runs);
    meter.measure([&](std::size_t i) { return pixels1[i] == pixels2[i]; });
  };

  BENCHMARK_ADVANCED("Pixel<N>: operator<")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
    auto rand_eng = get_prng();
    const auto pixels1 = random_sample_pixels(rand_eng, pixels_pop, num_runs);
    const auto pixels2 = random_sample_pixels(rand_eng, pixels_pop, num_runs);
    meter.measure([&](std::size_t i) { return pixels1[i] < pixels2[i]; });
  };

  BENCHMARK_ADVANCED("Pixel<N>: chrom()")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
    auto rand_eng = get_prng();
    const auto pixels = random_sample_pixels(rand_eng, pixels_pop, num_runs);
    meter.measure([&](std::size_t i) { return pixels[i].coords.bin1.chrom(); });
  };

  BENCHMARK_ADVANCED("Pixel<N>: chrom_id()")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
    auto rand_eng = get_prng();
    const auto pixels = random_sample_pixels(rand_eng, pixels_pop, num_runs);
    meter.measure([&](std::size_t i) { return pixels[i].coords.bin1.chrom().id(); });
  };

  BENCHMARK_ADVANCED("Pixel<N>: bin()")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = conditional_static_cast<std::size_t>(meter.runs());
    auto rand_eng = get_prng();
    const auto pixels = random_sample_pixels(rand_eng, pixels_pop, num_runs);
    meter.measure([&](std::size_t i) { return pixels[i].coords.bin1; });
  };
}

}  // namespace coolerpp::benchmark
