// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/chromosome.hpp"

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <string_view>
#include <vector>

namespace coolerpp::test::chromosome {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ChromosomeSet ctors", "[chromosome][short]") {
  // clang-format off
  const std::array<Chromosome, 3> expected{
      Chromosome{"chr1", 50001},
      Chromosome{"chr2", 25017},
      Chromosome{"chr3", 10000}
  };
  // clang-format on

  constexpr std::array<std::string_view, 3> expected_names{"chr1", "chr2", "chr3"};
  constexpr std::array<std::uint32_t, 3> expected_sizes{50001, 25017, 10000};

  SECTION("ctor w/ iterator of chromosomes") {
    const ChromosomeSet chroms(expected.begin(), expected.end());

    CHECK(chroms.size() == 3);
  }

  SECTION("ctor w/ iterator of chrom names and sizes") {
    const ChromosomeSet chroms(expected_names.begin(), expected_names.end(),
                               expected_sizes.begin());

    CHECK(chroms.size() == 3);
  }

  SECTION("ctor w/ iterator of chromosomes (duplicates)") {
    std::vector<Chromosome> expected_(expected.begin(), expected.end());
    expected_.push_back(expected.back());

    CHECK_THROWS_WITH(ChromosomeSet(expected_.begin(), expected_.end()),
                      Catch::Matchers::ContainsSubstring("found duplicate chromosome"));
  }

  SECTION("ctor w/ iterator of chrom names and sizes") {
    std::vector<std::string_view> expected_names_(expected_names.begin(), expected_names.end());
    std::vector<std::uint32_t> expected_sizes_(expected_sizes.begin(), expected_sizes.end());

    expected_names_.push_back(expected_names.back());
    expected_sizes_.push_back(expected_sizes.back());

    CHECK_THROWS_WITH(
        ChromosomeSet(expected_names_.begin(), expected_names_.end(), expected_sizes_.begin()),
        Catch::Matchers::ContainsSubstring("found duplicate chromosome"));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ChromosomeSet lookups", "[chromosome][short]") {
  // clang-format off
  const ChromosomeSet chroms{{
      Chromosome{"chr1", 50001},
      Chromosome{"chr2", 25017},
      Chromosome{"chr3", 10000}}
  };
  // clang-format on

  SECTION("contains") {
    CHECK(chroms.contains(Chromosome{"chr1", 50001}));
    CHECK(chroms.contains(0));
    CHECK(chroms.contains("chr1"));

    CHECK_FALSE(chroms.contains(Chromosome{"chr0", 50001}));
    CHECK_FALSE(chroms.contains(7));
    CHECK_FALSE(chroms.contains("chr0"));
    CHECK_FALSE(chroms.contains(""));
  }

  SECTION("at") {
    CHECK(chroms.at(0) == Chromosome{"chr1", 50001});
    CHECK(chroms.at("chr1") == Chromosome{"chr1", 50001});

    CHECK_THROWS_AS(chroms.at(3), std::out_of_range);
    CHECK_THROWS_AS(chroms.at("chr0"), std::out_of_range);
  }

  SECTION("get_id") {
    CHECK(chroms.get_id("chr1") == 0);
    CHECK(chroms.get_id(Chromosome{"chr3", 10000}) == 2);

    CHECK_THROWS_AS(chroms.get_id("a"), std::out_of_range);
  }
}
}  // namespace coolerpp::test::chromosome
