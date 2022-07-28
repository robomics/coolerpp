// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "coolerpp/coolerpp.hpp"
#include "coolerpp/test/self_deleting_folder.hpp"

namespace coolerpp::test {
inline const SelfDeletingFolder testdir{true};            // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace coolerpp::test

namespace coolerpp::test::coolerpp {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Coolerpp: format checking", "[cooler][short]") {
  SECTION("test .cool") {
    const auto path = datadir / "cooler_test_file.cool";
    CHECK(utils::is_cooler(path.string()));
    CHECK(!utils::is_multires_file(path.string()));
    CHECK(!utils::is_scool_file(path.string()));
  }

  SECTION("test .mcool") {
    const auto path = datadir / "multires_cooler_test_file.mcool";
    constexpr auto suffix{"::/resolutions/400000"};

    CHECK(!utils::is_cooler(path.string()));
    CHECK(utils::is_multires_file(path.string()));
    CHECK(!utils::is_scool_file(path.string()));
    CHECK(utils::is_cooler(path.string() + suffix));
  }

  SECTION("test .scool") {
    const auto path = datadir / "single_cell_cooler_test_file.scool";
    constexpr auto suffix{"::/cells/GSM2687248_41669_ACAGTG-R1-DpnII.100000.cool"};

    CHECK(!utils::is_cooler(path.string()));
    CHECK(!utils::is_multires_file(path.string()));
    CHECK(utils::is_scool_file(path.string()));
    CHECK(utils::is_cooler(path.string() + suffix));
  }

  SECTION("test with empty .h5 file") {
    const auto path = datadir / "empty_test_file.h5";
    CHECK(!utils::is_cooler(path.string()));
    CHECK(!utils::is_multires_file(path.string()));
    CHECK(!utils::is_scool_file(path.string()));
  }

  SECTION("test with nonexistent file") {
    const auto invalid_path = datadir / "void.nonexistent";
    CHECK_THROWS_WITH(utils::is_cooler(invalid_path.string()),
                      Catch::Matchers::ContainsSubstring("Unable to open file"));
    CHECK_THROWS_WITH(utils::is_multires_file(invalid_path.string()),
                      Catch::Matchers::ContainsSubstring("Unable to open file"));
    CHECK_THROWS_WITH(utils::is_scool_file(invalid_path.string()),
                      Catch::Matchers::ContainsSubstring("Unable to open file"));
  }

  SECTION("test corrupted .cool") {
    SECTION("missing format attribute") {
      const auto path = datadir / "invalid_coolers/missing_format_attr.cool";
      CHECK(utils::is_cooler(path.string()).missing_or_invalid_format_attr);
    }
    SECTION("invalid format attribute") {
      const auto path = datadir / "invalid_coolers/invalid_format_attr.cool";
      CHECK(utils::is_cooler(path.string()).missing_or_invalid_format_attr);
    }
  }

  SECTION("test corrupted .mcool") {
    // This file is missing group /resolutions/400000/pixels
    const auto path = datadir / "invalid_coolers/missing_pixels_group.mcool";
    const auto status = utils::is_multires_file(path.string());

    CHECK(!status);
    CHECK(status.is_hdf5);
    CHECK(!status.is_multires_file);
    CHECK(!status.missing_or_invalid_format_attr);
    CHECK(!status.missing_or_invalid_bin_type_attr);
    CHECK(status.uri == path.string());
    CHECK(status.missing_groups.empty());

    REQUIRE(status.invalid_resolutions.size() == 1);
    const auto& invalid_res = status.invalid_resolutions.front();

    const auto corrupted_uri_expected =
        fmt::format(FMT_STRING("{}::/resolutions/400000"), path.string());
    CHECK(invalid_res.uri == corrupted_uri_expected);
    CHECK(!invalid_res.is_cooler);
    REQUIRE(invalid_res.missing_groups.size() == 1);
    CHECK(invalid_res.missing_groups.front() == "pixels");
  }

  SECTION("test corrupted .scool") {
    // In this file, the number of groups under /cells and number of cells from ncells attribute
    // mismatch
    const auto path = datadir / "invalid_coolers/invalid_ncells_attribute.scool";
    const auto status = utils::is_scool_file(path.string());

    CHECK(!status);
    CHECK(status.is_hdf5);
    CHECK(!status.is_scool_file);
    CHECK(!status.missing_or_invalid_format_attr);
    CHECK(!status.missing_or_invalid_bin_type_attr);
    CHECK(status.uri == path.string());
    CHECK(status.missing_groups.empty());
    CHECK(status.unexpected_number_of_cells);
    CHECK(status.invalid_cells.empty());
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Coolerpp: file ctors", "[cooler][short]") {
  SECTION("open .cool") {
    const auto path = datadir / "cooler_test_file.cool";
    const auto f = File::open_read_only(path.string());

    CHECK(f.path() == path);
    CHECK(f.uri() == path);
    CHECK(f.bin_size() == 100'000);
    CHECK(f.chromosomes().size() == 20);
    CHECK(f.bins().size() == 26'398);
    CHECK(f.has_pixel_of_type<std::int32_t>());
  }

  SECTION("open .mcool") {
    const auto path = datadir / "single_cell_cooler_test_file.scool";

    SECTION("missing suffix") {
      CHECK_THROWS_WITH(
          File::open_read_only(path.string()),
          Catch::Matchers::ContainsSubstring("does not look like a valid Cooler file") &&
              Catch::Matchers::ContainsSubstring("missing_groups=[pixels, indexes]"));
    }

    SECTION("with suffix") {
      constexpr auto suffix{"::/cells/GSM2687248_41669_ACAGTG-R1-DpnII.100000.cool"};
      const auto f = File::open_read_only(path.string() + suffix);

      CHECK(f.path() == path);
      CHECK(f.uri() == path.string() + suffix);
    }
  }

  SECTION("open .scool") {
    const auto path = datadir / "multires_cooler_test_file.mcool";
    SECTION("missing suffix") {
      CHECK_THROWS_WITH(
          File::open_read_only(path.string()),
          Catch::Matchers::ContainsSubstring("does not look like a valid Cooler file") &&
              Catch::Matchers::ContainsSubstring("missing_groups=[chroms, bins, pixels, indexes]"));
    }

    SECTION("with suffix") {
      constexpr auto suffix = "::/resolutions/400000";

      const auto f = File::open_read_only(path.string() + suffix);
      CHECK(f.path() == path);
      CHECK(f.uri() == path.string() + suffix);
    }
  }

  SECTION("open empty .h5") {
    const auto path = datadir / "empty_test_file.h5";
    CHECK_THROWS_WITH(File::open_read_only(path.string()),
                      Catch::Matchers::ContainsSubstring("does not look like a valid Cooler file"));
  }

  SECTION("non existent") {
    const auto path = datadir / "cooler_test_file.cool.nonexistent";
    CHECK_THROWS_WITH(File::open_read_only(path.string()),
                      Catch::Matchers::ContainsSubstring("Unable to open file"));
  }

  SECTION("open corrupted .cool") {
    SECTION("corrupted bin table") {
      const auto path = datadir / "invalid_coolers/corrupted_bins.cool";
      CHECK_THROWS_WITH(File::open_read_only(path.string()),
                        Catch::Matchers::ContainsSubstring("Datasets have inconsistent sizes") &&
                            Catch::Matchers::ContainsSubstring("bins/chrom") &&
                            Catch::Matchers::ContainsSubstring("bins/start") &&
                            Catch::Matchers::ContainsSubstring("bins/end"));
    }

    SECTION("corrupted chrom table") {
      const auto path = datadir / "invalid_coolers/corrupted_chroms.cool";
      CHECK_THROWS_WITH(File::open_read_only(path.string()),
                        Catch::Matchers::ContainsSubstring("/chroms/name and") &&
                            Catch::Matchers::ContainsSubstring("/chroms/length shape mismatch"));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Coolerpp: init files", "[cooler][short]") {
  const ChromosomeSet chroms{Chromosome{"chr1", 10000}, Chromosome{"chr2", 5000}};

  SECTION(".cool") {
    const auto path = testdir() / "test_init.cool";
    constexpr std::uint32_t bin_size = 1000;
    std::ignore = File::create_new_cooler(path.string(), chroms, bin_size, true);
    CHECK(utils::is_cooler(path.string()));
  }

  SECTION(".mcool") {
    const auto path = testdir() / "test_init.mcool";
    constexpr std::array<std::uint32_t, 5> resolutions{10, 20, 30, 40, 50};
    init_mcool(path.string(), resolutions.begin(), resolutions.end(), true);

    for (const auto res : resolutions) {
      std::ignore = File::create_new_cooler(
          fmt::format(FMT_STRING("{}::/resolutions/{}"), path.string(), res), chroms, res);
    }

    CHECK(utils::is_multires_file(path.string()));
  }
  /*
    SECTION(".scool") {
      const auto path = (testdir() / "test_init.scool").string();
      constexpr std::array<std::string_view, 5> cell_ids{"1", "2", "3", "4", "5"};
      const std::array<std::pair<std::string_view, std::uint64_t>, 3> chroms{
          std::make_pair("chr1", 10000), std::make_pair("chr2", 5000), std::make_pair("chr3",
    1000)}; constexpr std::uint32_t bin_size = 50; init_scool(path, chroms.begin(), chroms.end(),
    cell_ids.begin(), cell_ids.end(), true);

      for (const auto id : cell_ids) {
        std::ignore =
            File::create_new_cooler(fmt::format(FMT_STRING("{}::/cells/{}"), path, id), bin_size);
      }

      CHECK(utils::is_scool_file(path));
    }
    */
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Coolerpp: read attributes", "[cooler][short]") {
  auto path = datadir / "cooler_test_file.cool";
  auto f = File::open_read_only(path.string());

  SECTION("bin size") { CHECK(f.bin_size() == 100'000); }

  SECTION("common attributes") {
    const auto& attrs = f.attributes();
    CHECK(attrs.bin_size == 100'000);
    CHECK(attrs.bin_type == "fixed");
    CHECK(attrs.creation_date == "2020-07-08T13:41:20.376258");
    CHECK(attrs.format == COOL_MAGIC);
    CHECK(attrs.format_url == "https://github.com/mirnylab/cooler");
    CHECK(attrs.format_version == 3);
    CHECK(attrs.generated_by == "cooler-0.8.8-dev");
    CHECK(attrs.assembly == "unknown");
    CHECK(attrs.metadata == "{}");
    CHECK(attrs.nbins == 26398);
    CHECK(attrs.nchroms == 20);
    CHECK(attrs.nnz == 107041);
    CHECK(attrs.storage_mode == "symmetric-upper");
    std::visit([](auto& sum) { CHECK(sum == 395465); }, attrs.sum);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Coolerpp: read/write chromosomes", "[cooler][short]") {
  const auto path = (testdir() / "test_write_chroms.cool").string();

  constexpr std::uint32_t bin_size = 5000;
  const ChromosomeSet chroms{Chromosome{"chr1", 50001}, Chromosome{"chr2", 25017},
                             Chromosome{"chr3", 10000}};

  {
    auto f = File::create_new_cooler(path, chroms, bin_size, true);
    CHECK(chroms == f.chromosomes());
  }

  const auto f = File::open_read_only(path, false);
  CHECK(chroms == f.chromosomes());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Coolerpp: read/write bin table", "[cooler][short]") {
  const auto path = (testdir() / "test_write_bin_table.cool").string();

  const ChromosomeSet chroms{Chromosome{"chr1", 50001}, Chromosome{"chr2", 25017},
                             Chromosome{"chr3", 10000}};

  constexpr std::uint32_t bin_size = 5000;
  const BinTableLazy table(chroms, bin_size);

  { auto f = File::create_new_cooler(path, chroms, bin_size, true); }

  auto f = File::open_read_only(path);

  auto start_it = f.dataset("bins/start").begin<std::uint32_t>();
  auto end_it = f.dataset("bins/end").begin<std::uint32_t>();

  REQUIRE(start_it != f.dataset("bins/start").end<std::uint32_t>());
  REQUIRE(end_it != f.dataset("bins/end").end<std::uint32_t>());

  for (const auto [_, start, end] : table) {
    CHECK(*start_it++ == start);
    CHECK(*end_it++ == end);
  }

  CHECK(start_it == f.dataset("bins/start").end<std::uint32_t>());
  CHECK(end_it == f.dataset("bins/end").end<std::uint32_t>());
}

}  // namespace coolerpp::test::coolerpp
