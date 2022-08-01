// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/pixel_selector.hpp"

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <random>

#include "coolerpp/coolerpp.hpp"
#include "coolerpp/test/self_deleting_folder.hpp"

namespace coolerpp::test {
inline const SelfDeletingFolder testdir{true};            // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace coolerpp::test

namespace coolerpp::test::pixel_selector {

template <class N>
static std::size_t generate_test_data(const std::filesystem::path& path,
                                      const ChromosomeSet& chroms, std::uint32_t bin_size) {
  auto f = File::create_new_cooler<N>(path.string(), chroms, bin_size, true);

  const auto num_bins = f.bins().size();

  std::vector<Pixel<N>> pixels;

  N n = 0;
  for (std::uint64_t i = 0; i < num_bins; ++i) {
    for (std::uint64_t j = i; j < num_bins; ++j) {
      pixels.emplace_back(Pixel<N>{PixelCoordinates{f.bins(), i, j}, n++});
    }
  }
  f.append_pixels(pixels.begin(), pixels.end());
  return pixels.size();
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Pixel selector: query", "[pixel_selector][short]") {
  const auto path = testdir() / "pixel_selector_devel.cool";

  const ChromosomeSet chroms{Chromosome{"chr1", 1000}};
  constexpr std::uint32_t bin_size = 10;
  using T = std::uint32_t;

  const auto expected_nnz = generate_test_data<T>(path, chroms, bin_size);

  auto f = File::open_read_only(path.string());
  const std::vector<Pixel<T>> expected_pixels(f.begin<T>(), f.end<T>());
  REQUIRE(expected_pixels.size() == expected_nnz);
  SECTION("query overlaps chrom start") {
    auto selector = f.fetch<T>("chr1:0-20");
    std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(std::distance(selector.begin(), selector.end()) == 3);
    CHECK(pixels[0].count == 0);
    CHECK(pixels[1].count == 1);
    CHECK(pixels[2].count == 100);
  }

  SECTION("query overlaps chrom end") {
    auto selector = f.fetch<T>("chr1:980-1000");
    std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(std::distance(selector.begin(), selector.end()) == 3);
    pixels = std::vector<Pixel<T>>(selector.begin(), selector.end());
    CHECK(pixels[0].count == 5047);
    CHECK(pixels[1].count == 5048);
    CHECK(pixels[2].count == 5049);
  }

  SECTION("query does not overlap chrom boundaries") {
    auto selector = f.fetch<T>("chr1:750-780");
    std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(std::distance(selector.begin(), selector.end()) == 6);
    pixels = std::vector<Pixel<T>>(selector.begin(), selector.end());
    CHECK(pixels[0].count == 4725);
    CHECK(pixels[1].count == 4726);
    CHECK(pixels[2].count == 4727);
    CHECK(pixels[3].count == 4750);
    CHECK(pixels[4].count == 4751);
    CHECK(pixels[5].count == 4774);
  }

  SECTION("query does not line up with bins") {
    auto selector = f.fetch<T>("chr1:901-927");
    std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(std::distance(selector.begin(), selector.end()) == 6);
    pixels = std::vector<Pixel<T>>(selector.begin(), selector.end());
    CHECK(pixels[0].count == 4995);
    CHECK(pixels[1].count == 4996);
    CHECK(pixels[2].count == 4997);
    CHECK(pixels[3].count == 5005);
    CHECK(pixels[4].count == 5006);
    CHECK(pixels[5].count == 5014);
  }

  SECTION("large query") {
    auto selector = f.fetch<T>("chr1:75-975");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 4186);

    const auto sum = std::accumulate(
        selector.begin(), selector.end(), T(0),
        [&](T accumulator, const Pixel<T> pixel) { return accumulator + pixel.count; });

    CHECK(sum == 11852659);
  }

  SECTION("empty query") {
    auto selector = f.fetch<T>("chr1:0-0");
    CHECK(selector.begin() == selector.end());

    selector = f.fetch<T>("chr1:100-100");
    CHECK(selector.begin() == selector.end());

    selector = f.fetch<T>("chr1:1000-1000");
    CHECK(selector.begin() == selector.end());
  }

  SECTION("query spans 1 bin") {
    auto selector = f.fetch<T>("chr1:0-9");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 0);

    selector = f.fetch<T>("chr1:5-7");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 0);

    selector = f.fetch<T>("chr1:991-1000");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 5049);
  }

  SECTION("query spans 1bp") {
    auto selector = f.fetch<T>("chr1:0-1");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 0);

    selector = f.fetch<T>("chr1:12-13");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 100);

    selector = f.fetch<T>("chr1:999-1000");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 5049);
  }

  SECTION("query spans entire chromosome") {
    auto selector = f.fetch<T>("chr1");

    CHECK(std::distance(selector.begin(), selector.end()) == 5050);
    const auto sum = std::accumulate(
        selector.begin(), selector.end(), T(0),
        [&](T accumulator, const Pixel<T> pixel) { return accumulator + pixel.count; });
    CHECK(sum == 12748725);
  }

  SECTION("invalid queries") {
    CHECK_THROWS_WITH(f.fetch<T>(""), Catch::Matchers::Equals("query \"\" is malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr2:0-1"),
                      Catch::Matchers::ContainsSubstring("invalid chromosome"));

    CHECK_THROWS_WITH(f.fetch<T>(":0-1"), Catch::Matchers::ContainsSubstring("invalid chromosome"));
    CHECK_THROWS_WITH(f.fetch<T>("-:0-1"),
                      Catch::Matchers::ContainsSubstring("invalid chromosome"));
    CHECK_THROWS_WITH(f.fetch<T>("::0-1"),
                      Catch::Matchers::ContainsSubstring("invalid chromosome"));

    CHECK_THROWS_WITH(f.fetch<T>("chr1:"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1-"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:-"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1-0-1"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:0:1"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:01"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:-01"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:-01"), Catch::Matchers::ContainsSubstring("malformed"));

    CHECK_THROWS_WITH(f.fetch<T>("chr1:-1"),
                      Catch::Matchers::ContainsSubstring("missing start position"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:0-"),
                      Catch::Matchers::ContainsSubstring("missing end position"));

    CHECK_THROWS_WITH(f.fetch<T>("chr1:4294967296-0"),
                      Catch::Matchers::ContainsSubstring("invalid start position"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:0-4294967296"),
                      Catch::Matchers::ContainsSubstring("invalid end position"));

    CHECK_THROWS_WITH(f.fetch<T>("chr1:10-5"), Catch::Matchers::ContainsSubstring(
                                                   "start position is greater than end position"));
  }
}

}  // namespace coolerpp::test::pixel_selector
