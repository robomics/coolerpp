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
static std::ptrdiff_t generate_test_data(const std::filesystem::path& path,
                                         const ChromosomeSet& chroms, std::uint32_t bin_size) {
  auto f = File::create_new_cooler<N>(path.string(), chroms, bin_size, true);

  const auto num_bins = f.bins().size();

  std::vector<Pixel<N>> pixels;

  N n = 0;
  for (std::uint64_t i = 0; i < num_bins; ++i) {
    for (std::uint64_t j = i; j < num_bins; ++j) {
      pixels.emplace_back(Pixel<N>{PixelCoordinates{f.bins_ptr(), i, j}, n++});
    }
  }
  f.append_pixels(pixels.begin(), pixels.end());
  return static_cast<std::ptrdiff_t>(pixels.size());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Pixel selector: 1D queries", "[pixel_selector][short]") {
  const auto path1 = testdir() / "pixel_selector_devel.cool";

  const ChromosomeSet chroms{Chromosome{"chr1", 1000}};
  constexpr std::uint32_t bin_size = 10;
  using T = std::uint32_t;

  const auto expected_nnz = generate_test_data<T>(path1, chroms, bin_size);

  auto f = File::open_read_only(path1.string());
  REQUIRE(std::distance(f.begin<T>(), f.end<T>()) == expected_nnz);
  SECTION("query overlaps chrom start") {
    auto selector = f.fetch<T>("chr1:0-20");
    const std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(pixels.size() == 3);

    CHECK(pixels[0].count == 0);
    CHECK(pixels[1].count == 1);
    CHECK(pixels[2].count == 100);
  }

  SECTION("query overlaps chrom end") {
    auto selector = f.fetch<T>("chr1:980-1000");
    const std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(pixels.size() == 3);

    CHECK(pixels[0].count == 5047);
    CHECK(pixels[1].count == 5048);
    CHECK(pixels[2].count == 5049);
  }

  SECTION("query does not overlap chrom boundaries") {
    auto selector = f.fetch<T>("chr1:750-780");
    const std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(pixels.size() == 6);

    CHECK(pixels[0].count == 4725);
    CHECK(pixels[1].count == 4726);
    CHECK(pixels[2].count == 4727);
    CHECK(pixels[3].count == 4750);
    CHECK(pixels[4].count == 4751);
    CHECK(pixels[5].count == 4774);
  }

  SECTION("query does not line up with bins") {
    auto selector = f.fetch<T>("chr1:901-927");
    const std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(pixels.size() == 6);

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

    CHECK_THROWS_WITH(f.fetch<T>("chr1:0-0"),
                      Catch::Matchers::ContainsSubstring(
                          "end position should be greater than the start position"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:10-5"),
                      Catch::Matchers::ContainsSubstring(
                          "end position should be greater than the start position"));
  }
}
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Pixel selector: 2D queries", "[pixel_selector][short]") {
  using T = std::uint32_t;
  const auto path = datadir / "cooler_test_file.cool";
  auto f = File::open_read_only(path.string());

  SECTION("cis") {
    SECTION("valid") {
      auto selector = f.fetch<T>("1:5000000-5500000", "1:5000000-6500000");
      const std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
      REQUIRE(pixels.size() == 8);

      CHECK(pixels[0].count == 20);
      CHECK(pixels[1].count == 1);
      CHECK(pixels[2].count == 18);
      CHECK(pixels[3].count == 8);
      CHECK(pixels[4].count == 1);
      CHECK(pixels[5].count == 9);
      CHECK(pixels[6].count == 6);
      CHECK(pixels[7].count == 2);
    }

    SECTION("empty") {
      auto selector = f.fetch<T>("1:0-100000");
      CHECK(selector.begin() == selector.end());
    }
  }

  SECTION("trans") {
    SECTION("valid") {
      auto selector = f.fetch<T>("1:48000000-50000000", "4:30000000-35000000");
      const std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
      REQUIRE(pixels.size() == 6);

      CHECK(pixels[0].count == 1);
      CHECK(pixels[1].count == 3);
      CHECK(pixels[2].count == 1);
      CHECK(pixels[3].count == 3);
      CHECK(pixels[4].count == 7);
      CHECK(pixels[5].count == 1);
    }

    SECTION("empty") {
      auto selector = f.fetch<T>("1:0-50000", "2:0-50000");
      CHECK(selector.begin() == selector.end());
    }
  }
}

}  // namespace coolerpp::test::pixel_selector
