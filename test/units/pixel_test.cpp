// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/pixel.hpp"

#include <catch2/catch_test_macros.hpp>

#include "coolerpp/chromosome.hpp"

namespace coolerpp {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Pixel", "[pixel][short]") {
  const ChromosomeSet chroms{Chromosome{"chr1", 248'956'422}, Chromosome{"chr2", 242'193'529},
                             Chromosome{"chr3", 198'295'559}, Chromosome{"chr4", 190'214'555},
                             Chromosome{"chr5", 181'538'259}, Chromosome{"chr6", 170'805'979},
                             Chromosome{"chr9", 138'394'717}, Chromosome{"chr11", 135'086'622},
                             Chromosome{"chr12", 133'275'309}};
  constexpr std::uint32_t bin_size = 1;
  const auto bins = std::make_shared<const BinTableLazy>(chroms, bin_size);

  auto P = [&](std::string_view chrom1, std::string_view chrom2, std::uint32_t pos1,
               std::uint32_t pos2, std::uint32_t count = 0) {
    return Pixel<std::uint32_t>{PixelCoordinates{bins, chrom1, chrom2, pos1, pos2}, count};
  };

  SECTION("operator bool") {
    CHECK(!PixelCoordinates{});
    CHECK(!!P("chr1", "chr1", 0, 10));
  }

  SECTION("(dis)equality") {
    CHECK(P("chr1", "chr1", 0, 10) == P("chr1", "chr1", 0, 10));

    CHECK(P("chr1", "chr1", 0, 10) != P("chr1", "chr2", 0, 10));
    CHECK(P("chr1", "chr1", 0, 10) != P("chr2", "chr1", 0, 10));

    CHECK(P("chr1", "chr1", 0, 10) != P("chr1", "chr1", 1, 10));
    CHECK(P("chr1", "chr1", 0, 10) != P("chr1", "chr1", 0, 11));
  }

  SECTION("ordering") {
    CHECK(P("chr1", "chr1", 0, 0) < P("chr2", "chr2", 0, 0));
    CHECK(P("chr1", "chr1", 0, 0) <= P("chr2", "chr2", 0, 0));

    CHECK(P("chr1", "chr1", 0, 0) < P("chr1", "chr2", 0, 0));
    CHECK(P("chr1", "chr1", 0, 0) <= P("chr1", "chr2", 0, 0));

    CHECK(P("chr2", "chr2", 0, 0) > P("chr1", "chr1", 0, 0));
    CHECK(P("chr2", "chr2", 0, 0) >= P("chr1", "chr1", 0, 0));

    CHECK(P("chr1", "chr2", 0, 0) > P("chr1", "chr1", 0, 0));
    CHECK(P("chr1", "chr2", 0, 0) >= P("chr1", "chr1", 0, 0));

    CHECK(P("chr1", "chr1", 0, 0) < P("chr1", "chr1", 0, 1));
    CHECK(P("chr1", "chr1", 0, 0) < P("chr1", "chr1", 1, 0));
    CHECK(P("chr1", "chr1", 0, 0) <= P("chr1", "chr1", 0, 1));
    CHECK(P("chr1", "chr1", 0, 0) <= P("chr1", "chr1", 1, 0));

    CHECK(P("chr1", "chr1", 0, 1) > P("chr1", "chr1", 0, 0));
    CHECK(P("chr1", "chr1", 1, 0) > P("chr1", "chr1", 0, 0));
    CHECK(P("chr1", "chr1", 0, 1) >= P("chr1", "chr1", 0, 0));
    CHECK(P("chr1", "chr1", 1, 0) >= P("chr1", "chr1", 0, 0));
  }

  SECTION("sorting") {
    const std::vector<PixelCoordinates> coords{
        // clang-format off
        {bins, "chr1", "chr1",  10'000, 180'000},
        {bins, "chr1", "chr1",  10'000, 202'890'000},
        {bins, "chr1", "chr2",  10'000, 113'590'000},
        {bins, "chr1", "chr4",  10'000, 52'880'000},
        {bins, "chr1", "chr5",  10'000, 230'000},
        {bins, "chr1", "chr6",  10'000, 33'820'000},
        {bins, "chr1", "chr6",  10'000, 149'280'000},
        {bins, "chr1", "chr9",  10'000, 10'000},
        {bins, "chr1", "chr9",  10'000, 122'380'000},
        {bins, "chr1", "chr11", 40'000, 11'630'000},
        {bins, "chr1", "chr11", 40'000, 120'770'000},
        {bins, "chr1", "chr12", 40'000, 7'060'000},
        {bins, "chr1", "chr12", 40'000, 119'750'000},
        {bins, "chr2", "chr2",  10'000, 10'000},
        {bins, "chr2", "chr2",  10'000, 20'000},
        {bins, "chr2", "chr3",  10'000, 99'320'000},
        {bins, "chr2", "chr3",  10'000, 101'660'000},
        // clang-format on
    };

    CHECK(std::is_sorted(coords.begin(), coords.end()));
  }
}

}  // namespace coolerpp
