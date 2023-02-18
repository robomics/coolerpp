// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/index.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "coolerpp/bin_table.hpp"
#include "coolerpp/chromosome.hpp"

namespace coolerpp::test {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}

namespace coolerpp::test::index {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Index: ctor", "[index][short]") {
  constexpr std::uint32_t bin_size = 100;
  const auto bins = std::make_shared<const BinTableLazy>(
      ChromosomeSet{Chromosome{"chr1", 10001}, Chromosome{"chr2", 5000}}, bin_size);

  const Index idx(bins);

  CHECK(idx.bin_size() == bin_size);
  CHECK(idx.num_chromosomes() == 2);
  CHECK(idx.size() == 151);

  CHECK(idx.size("chr1") == 101);
  CHECK(idx.size(0) == 101);

  CHECK(idx.size("chr2") == 50);
  CHECK(idx.size(1) == 50);

  CHECK_THROWS_WITH(idx.size("chr3"), Catch::Matchers::Equals("chromosome \"chr3\" not found"));
  CHECK_THROWS_WITH(idx.size(99), Catch::Matchers::Equals("chromosome with id 99 not found"));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Index: offset setters and getters", "[index][short]") {
  constexpr std::uint32_t bin_size = 10;
  const auto bins =
      std::make_shared<const BinTableLazy>(ChromosomeSet{Chromosome{"chr1", 100}}, bin_size);

  constexpr auto fill_value = std::numeric_limits<std::uint64_t>::max();

  Index idx(bins);

  SECTION("by pos") {
    idx.set_offset_by_pos("chr1", 22, 1);
    idx.set_offset_by_pos(0, 55, 1);

    for (std::uint32_t pos = 0; pos < 100; ++pos) {
      const auto row_idx = conditional_static_cast<std::size_t>(pos / bin_size);
      const std::uint64_t expected = (row_idx == 2 || row_idx == 5) ? 1 : fill_value;

      CHECK(idx.get_offset_by_row_idx(0, row_idx) == expected);

      CHECK(idx.get_offset_by_pos("chr1", pos) == expected);
      CHECK(idx.get_offset_by_pos(0, pos) == expected);
    }
  }

  SECTION("by row idx") {
    idx.set_offset_by_row_idx(0, 2, 1);
    idx.set_offset_by_row_idx(0, 5, 1);

    for (std::uint32_t pos = 0; pos < 100; ++pos) {
      const auto row_idx = conditional_static_cast<std::size_t>(pos / bin_size);
      const std::size_t expected = (row_idx == 2 || row_idx == 5) ? 1 : fill_value;

      CHECK(idx.get_offset_by_row_idx(0, row_idx) == expected);

      CHECK(idx.get_offset_by_pos("chr1", pos) == expected);
      CHECK(idx.get_offset_by_pos(0, pos) == expected);
    }
  }

  SECTION("by bin ID") {
    idx.set_offset_by_bin_id(9, 9);
    CHECK(idx.get_offset_by_pos("chr1", 99) == 9);
    CHECK(idx.get_offset_by_bin_id(9) == 9);
  }

  SECTION("out of bound access") {
    CHECK_THROWS_WITH(idx.get_offset_by_pos("chr1", 999),
                      Catch::Matchers::ContainsSubstring("row maps outside of chromosome"));

    CHECK_THROWS_WITH(idx.get_offset_by_row_idx(0, 999),
                      Catch::Matchers::ContainsSubstring("row maps outside of chromosome"));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Index: iterator", "[index][short]") {
  constexpr std::uint32_t bin_size = 1000;
  const auto bins = std::make_shared<const BinTableLazy>(
      ChromosomeSet{Chromosome{"chr1", 10001}, Chromosome{"chr2", 5000}}, bin_size);

  // Assume there are 10 pixels per row
  constexpr std::array<std::size_t, 11> chr1_offsets{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  constexpr std::array<std::size_t, 5> chr2_offsets{110, 120, 130, 140, 150};

  Index idx(bins);
  for (std::size_t i = 0; i < chr1_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(0, i, chr1_offsets[i]);
  }

  for (std::size_t i = 0; i < chr2_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(1, i, chr2_offsets[i]);
  }

  auto first_offset = idx.begin();
  const auto last_offset = idx.end();
  REQUIRE(first_offset != last_offset);

  for (std::size_t i = 0; i < idx.size(); ++i) {
    CHECK(*first_offset++ == conditional_static_cast<std::uint64_t>(i * 10));
  }

  REQUIRE(first_offset != last_offset);
  CHECK(*first_offset == 0);
  idx.finalize(16);
  CHECK(*first_offset++ == idx.size());
  CHECK(first_offset == last_offset);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Index: validation", "[index][short]") {
  constexpr std::uint32_t bin_size = 1000;
  const auto bins = std::make_shared<const BinTableLazy>(
      ChromosomeSet{Chromosome{"chr1", 10001}, Chromosome{"chr2", 5000}}, bin_size);

  // Assume there are 10 pixels per row
  constexpr std::array<std::size_t, 11> chr1_offsets{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  constexpr std::array<std::size_t, 5> chr2_offsets{110, 120, 130, 140, 150};

  Index idx(bins);
  for (std::size_t i = 0; i < chr1_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(0, i, chr1_offsets[i]);
  }

  for (std::size_t i = 0; i < chr2_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(1, i, chr2_offsets[i]);
  }

  SECTION("valid index") { CHECK_NOTHROW(idx.validate()); }
  SECTION("first offset is not zero") {
    idx.set_offset_by_row_idx(0, 0, 1);
    CHECK_THROWS_WITH(idx.validate(),
                      Catch::Matchers::ContainsSubstring("first offset is not zero"));
  }
  SECTION("offsets for adjacent chromosomes are not in ascending order") {
    idx.set_offset_by_row_idx(1, 0, 99);
    CHECK_THROWS_WITH(idx.validate(), Catch::Matchers::ContainsSubstring(
                                          "offset for bin chr2:0-1000 should be >= 100, found 99"));
  }

  SECTION("offsets are not sorted") {
    idx.set_offset_by_row_idx(1, 2, 150);
    CHECK_THROWS_WITH(idx.validate(),
                      Catch::Matchers::ContainsSubstring("offsets are not in ascending order"));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Index: compute chromosome offsets", "[index][short]") {
  constexpr std::uint32_t bin_size = 1000;
  const auto bins = std::make_shared<const BinTableLazy>(
      ChromosomeSet{Chromosome{"chr1", 10001}, Chromosome{"chr2", 5000}}, bin_size);

  // Assume there are 10 pixels per row
  constexpr std::array<std::size_t, 11> chr1_offsets{10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  constexpr std::array<std::size_t, 5> chr2_offsets{110, 120, 130, 140, 150};

  Index idx(bins);

  for (std::size_t i = 0; i < chr1_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(0, i, chr1_offsets[i]);
  }

  for (std::size_t i = 0; i < chr2_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(1, i, chr2_offsets[i]);
  }

  const auto chrom_offsets = idx.compute_chrom_offsets();
  REQUIRE(chrom_offsets.size() == bins->num_chromosomes() + 1);

  CHECK(chrom_offsets[0] == 0);
  CHECK(chrom_offsets[1] == chr1_offsets.size());
  CHECK(chrom_offsets[2] == chr1_offsets.size() + chr2_offsets.size());
}

}  // namespace coolerpp::test::index
