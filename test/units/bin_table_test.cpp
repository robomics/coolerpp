// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/bin_table.hpp"

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace coolerpp::test::bin_table {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BinTable", "[bin-table][short]") {
  constexpr std::uint32_t bin_size = 5000;
  // clang-format off
  const BinTable table({
      Chromosome{0, "chr1", 50001},
      Chromosome{1, "chr2", 25017},
      Chromosome{2, "chr3", 10000}},
      bin_size);
  // clang-format on

  SECTION("stats") {
    CHECK(table.size() == 11 + 6 + 2);
    CHECK(table.num_chromosomes() == 3);
    CHECK(table.bin_size() == bin_size);
  }

  SECTION("at") {
    const BinTable expected{{Chromosome{1, "chr2", 25017}}, bin_size};

    CHECK(table.subset(Chromosome{1, "chr2", 25017}) == expected);
    CHECK(table.subset("chr2") == expected);
    CHECK(table.subset(1) == expected);
    CHECK(table.subset("chr1") != expected);

    CHECK_THROWS_AS(table.subset(Chromosome{4, "chr5", 1}), std::out_of_range);
    CHECK_THROWS_AS(table.subset("a"), std::out_of_range);
    CHECK_THROWS_AS(table.subset(10), std::out_of_range);
  }

  SECTION("bin id to coord") {
    const auto& chr1 = table.chromosomes().at("chr1");
    const auto& chr2 = table.chromosomes().at("chr2");

    CHECK(table.at(0) == Bin{chr1, 0, bin_size});
    CHECK(table.at(10) == Bin{chr1, 50000, 50001});

    CHECK(table.at(11) == Bin{chr2, 0, bin_size});

    CHECK_THROWS_AS(table.at(table.size()), std::out_of_range);
  }

  SECTION("coord to bin id") {
    const auto& chr2 = table.chromosomes().at("chr2");

    CHECK(table.map_to_bin_id(0, 7500) == 1);
    CHECK(table.map_to_bin_id("chr1", 50000) == 10);
    CHECK(table.map_to_bin_id(chr2, 10) == 11);

    CHECK_THROWS_AS(table.map_to_bin_id("a", 0), std::out_of_range);
    CHECK_THROWS_AS(table.map_to_bin_id("chr1", 99999), std::out_of_range);
    CHECK_THROWS_AS(table.map_to_bin_id(chr2, 99999), std::out_of_range);
    CHECK_THROWS_AS(table.map_to_bin_id(1, 99999), std::out_of_range);
  }

  SECTION("iterators") {
    const auto& chr1 = table.chromosomes().at("chr1");
    const auto& chr2 = table.chromosomes().at("chr2");
    const auto& chr3 = table.chromosomes().at("chr3");

    // clang-format off
    const std::array<Bin, 19> expected{
       Bin{chr1, 0, 5000},
       Bin{chr1, 5000, 10000},
       Bin{chr1, 10000, 15000},
       Bin{chr1, 15000, 20000},
       Bin{chr1, 20000, 25000},
       Bin{chr1, 25000, 30000},
       Bin{chr1, 30000, 35000},
       Bin{chr1, 35000, 40000},
       Bin{chr1, 40000, 45000},
       Bin{chr1, 45000, 50000},
       Bin{chr1, 50000, 50001},
       Bin{chr2, 0, 5000},
       Bin{chr2, 5000, 10000},
       Bin{chr2, 10000, 15000},
       Bin{chr2, 15000, 20000},
       Bin{chr2, 20000, 25000},
       Bin{chr2, 25000, 25017},
       Bin{chr3, 0, 5000},
       Bin{chr3, 5000, 10000}
    };
    // clang-format on

    REQUIRE(table.size() == expected.size());

    SECTION("forward") {
      auto first_bin = table.begin();
      auto last_bin = table.end();

      // NOLINTNEXTLINE
      for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(*first_bin++ == expected[i]);
      }

      CHECK(first_bin == last_bin);
    }

    SECTION("backward") {
      auto first_bin = table.begin();
      auto last_bin = table.end();

      // NOLINTNEXTLINE
      for (std::size_t i = expected.size(); i != 0; --i) {
        CHECK(*(--last_bin) == expected[i - 1]);
      }

      CHECK(first_bin == last_bin);
    }
  }

  SECTION("concretize") {
    const auto concrete_table = table.concretize();

    REQUIRE(concrete_table.chroms.size() == table.size());

    std::size_t i = 0;
    for (const auto& bin : table) {
      CHECK(*concrete_table.chroms[i] == bin.chrom());
      CHECK(concrete_table.bin_starts[i] == bin.start());
      CHECK(concrete_table.bin_ends[i++] == bin.end());
    }
  }
}
}  // namespace coolerpp::test::bin_table
