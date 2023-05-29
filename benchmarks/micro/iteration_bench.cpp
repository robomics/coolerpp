// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <array>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <numeric>

#include "coolerpp/benchmark/common.hpp"
#include "coolerpp/coolerpp.hpp"
#include "coolerpp/pixel.hpp"

namespace coolerpp::benchmark {

inline const std::filesystem::path datadir{"benchmarks/data"};  // NOLINT(cert-err58-cpp)

TEST_CASE("Dataset iteration: benchmark", "[iteration][bench]") {
  const auto test_file = datadir / "4DNFI9FVHJZQ.0.9.1.mcool::/resolutions/1000000";

  BENCHMARK_ADVANCED("Dataset iteration: operator++")
  (Catch::Benchmark::Chronometer meter) {
    const auto clr = File::open_read_only(test_file.string());

    auto dset = clr.dataset("pixels/count");

    const std::size_t iters = 1'000'000;
    REQUIRE(dset.size() > iters);  // NOLINT

    meter.measure([&]() {
      std::int64_t n{};
      auto first = dset.begin<std::int64_t>();
      for (std::size_t i = 0; i < iters; ++i) {
        n += *++first;
      }
      return n;
    });
  };

  BENCHMARK_ADVANCED("Dataset iteration: operator++(int)")
  (Catch::Benchmark::Chronometer meter) {
    const auto clr = File::open_read_only(test_file.string());

    auto dset = clr.dataset("pixels/count");

    const std::size_t iters = 1'000'000;
    REQUIRE(dset.size() > iters);  // NOLINT

    meter.measure([&]() {
      std::int64_t n{};
      auto first = dset.begin<std::int64_t>();
      for (std::size_t i = 0; i < iters; ++i) {
        n += *first++;
      }
      return n;
    });
  };
}

TEST_CASE("Pixel<N> iteration: benchmark", "[iteration][bench]") {
  const auto test_file = datadir / "4DNFI9FVHJZQ.0.9.1.mcool::/resolutions/1000000";

  BENCHMARK_ADVANCED("Pixel<N> iteration: operator++")
  (Catch::Benchmark::Chronometer meter) {
    const auto clr = File::open_read_only(test_file.string());

    const std::int64_t iters = 1'000'000;
    REQUIRE(*clr.attributes().nnz > iters);  // NOLINT

    meter.measure([&]() {
      std::int64_t n{};
      auto first = clr.begin<std::int32_t>();
      for (std::int64_t i = 0; i < iters; ++i) {
        n += (*++first).count;
      }
      return n;
    });
  };

  BENCHMARK_ADVANCED("Pixel<N> iteration: operator++(int)")
  (Catch::Benchmark::Chronometer meter) {
    const auto clr = File::open_read_only(test_file.string());

    const std::int64_t iters = 1'000'000;
    REQUIRE(*clr.attributes().nnz > iters);  // NOLINT

    meter.measure([&]() {
      std::int64_t n{};
      auto first = clr.begin<std::int32_t>();
      for (std::int64_t i = 0; i < iters; ++i) {
        n += (*first++).count;
      }
      return n;
    });
  };
}

template <std::size_t chunk_size, std::size_t iters>
static void pixel_iterator_chunk_size_bench() {
  BENCHMARK_ADVANCED("File::iterator CHUNK_SIZE=" + std::to_string(chunk_size))
  (Catch::Benchmark::Chronometer meter) {
    const auto test_file = datadir / "4DNFI9FVHJZQ.0.9.1.mcool::/resolutions/100000";
    const auto clr = File::open_read_only(test_file.string());

    REQUIRE(*clr.attributes().nnz > static_cast<std::int64_t>(iters));  // NOLINT

    meter.measure([&]() {
      std::int64_t n{};
      auto first = clr.begin<std::int32_t, chunk_size>();
      for (std::size_t i = 0; i < iters; ++i) {
        n += (*++first).count;
      }
      return n;
    });
  };
}

TEST_CASE("File::iterator chunk_size: benchmark", "[iteration][bench]") {
  constexpr std::size_t MAX_ITERS = 1024 * 1024;
  pixel_iterator_chunk_size_bench<256, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<512, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<1024, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<2 * 1024, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<4 * 1024, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<8 * 1024, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<16 * 1024, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<32 * 1024, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<64 * 1024, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<128 * 1024, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<256 * 1024, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<512 * 1024, MAX_ITERS>();
  pixel_iterator_chunk_size_bench<1024 * 1024, MAX_ITERS>();
}

template <std::size_t chunk_size>
static void pixel_selector_iterator_chunk_size_bench(std::string query1, std::string query2) {
  BENCHMARK_ADVANCED(fmt::format(FMT_STRING("PixelSelector<N>::iterator {}:{} CHUNK_SIZE={}"),
                                 query1, query2, chunk_size))
  (Catch::Benchmark::Chronometer meter) {
    const auto test_file = datadir / "4DNFI9FVHJZQ.0.9.1.mcool::/resolutions/100000";
    const auto clr = File::open_read_only(test_file.string());

    // Returns ~ 4M pixels
    const auto sel = clr.fetch<std::int32_t>(query1, query2);

    meter.measure([&]() {
      return std::accumulate(
          sel.begin(), sel.end(), std::int64_t(0),
          [](auto accumulator, const auto& pixel) { return accumulator + pixel.count; });
    });
  };
}

TEST_CASE("PixelSelector<N>::iterator (CIS large): benchmark", "[iteration][bench]") {
  const std::string query = "chr1:0-50000000";  // returns ~120k pixels
  pixel_selector_iterator_chunk_size_bench<2 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<4 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<8 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<16 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<32 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<64 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<128 * 1024>(query, query);
}

TEST_CASE("PixelSelector<N>::iterator (CIS small): benchmark", "[iteration][bench]") {
  const std::string query = "chr1:10000000-15000000";  // returns ~1k pixels
  pixel_selector_iterator_chunk_size_bench<256>(query, query);
  pixel_selector_iterator_chunk_size_bench<512>(query, query);
  pixel_selector_iterator_chunk_size_bench<1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<2 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<4 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<8 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<16 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<32 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<64 * 1024>(query, query);
  pixel_selector_iterator_chunk_size_bench<128 * 1024>(query, query);
}

TEST_CASE("PixelSelector<N>::iterator (TRANS large): benchmark", "[iteration][bench]") {
  // returns ~300k pixels
  const std::string query1 = "chr1:0-50000000";
  const std::string query2 = "chr2:0-50000000";
  pixel_selector_iterator_chunk_size_bench<2 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<4 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<8 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<16 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<32 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<64 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<128 * 1024>(query1, query2);
}

TEST_CASE("PixelSelector<N>::iterator (TRANS small): benchmark", "[iteration][bench]") {
  // returns ~2k pixels
  const std::string query1 = "chr1:10000000-15000000";
  const std::string query2 = "chr2:50000000-55000000";
  pixel_selector_iterator_chunk_size_bench<256>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<512>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<2 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<4 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<8 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<16 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<32 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<64 * 1024>(query1, query2);
  pixel_selector_iterator_chunk_size_bench<128 * 1024>(query1, query2);
}

}  // namespace coolerpp::benchmark
