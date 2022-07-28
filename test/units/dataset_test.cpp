// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/dataset.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <random>
#include <set>

#include "coolerpp/group.hpp"
#include "coolerpp/test/self_deleting_folder.hpp"

namespace coolerpp::test {
inline const SelfDeletingFolder testdir{true};            // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace coolerpp::test

namespace coolerpp::test::dataset {
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Dataset - read", "[dataset][short]") {
  const auto path = datadir / "cooler_test_file.cool";
  RootGroup grp{HighFive::File(path.string()).getGroup("/")};

  SECTION("fixed str") {
    SECTION("vector") {
      constexpr std::array<std::string_view, 3> expected{"1", "2", "3"};
      std::vector<std::string> buff{};

      Dataset{grp, "chroms/name"}.read(buff, expected.size());
      REQUIRE(buff.size() == 3);

      for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }
    }

    SECTION("atomic") {
      Dataset dset{grp, "chroms/name"};
      std::string buff;
      dset.read(buff, 9);
      CHECK(buff == "10");
      CHECK(dset.read_last<std::string>() == "X");
      CHECK(std::get<std::string>(dset.read_last()) == "X");
    }
  }

  SECTION("numeric") {
    using T = std::int32_t;
    constexpr std::array<T, 10> expected{0,       100'000, 200'000, 300'000, 400'000,
                                         500'000, 600'000, 700'000, 800'000, 900'000};

    constexpr std::size_t nnz_expected = 107'041;
    constexpr std::int32_t sum_expected = 395'465;

    SECTION("vector<T>") {
      std::vector<T> buff{};
      std::ignore = Dataset{grp, "bins/start"}.read(buff, expected.size());

      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }

      Dataset{grp, "pixels/count"}.read_all(buff);
      CHECK(buff.size() == nnz_expected);
      CHECK(std::accumulate(buff.begin(), buff.end(), 0) == sum_expected);
    }

    SECTION("variant buff") {
      internal::VariantBuffer vbuff{std::size_t(0), 0.0};
      std::ignore = Dataset{grp, "bins/start"}.read(vbuff, expected.size());
      const auto& buff = vbuff.get<T>();
      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }

      Dataset{grp, "pixels/count"}.read_all(vbuff);
      CHECK(vbuff.size<T>() == nnz_expected);
      CHECK(std::accumulate(vbuff.begin<T>(), vbuff.end<T>(), 0) == sum_expected);
    }

    SECTION("atomic") {
      Dataset dset{grp, "chroms/length"};
      std::uint64_t buff;
      dset.read(buff, 2);
      CHECK(buff == 159'599'783);

      CHECK(dset.read_last<std::int32_t>() == 166'650'296);
      CHECK(std::get<std::int32_t>(dset.read_last()) == 166'650'296);
    }

    SECTION("enum") { CHECK(Dataset{grp, "bins/chrom"}.read<std::uint32_t>(0) == 0); }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Dataset - write", "[dataset][short]") {
  const auto path = testdir() / "test_dataset_write.cool";
  RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};

  SECTION("fixed str") {
    SECTION("vector") {
      using BuffT = std::vector<std::string>;
      const BuffT expected{"s1", "this_is_a_relatively_long_string"};

      Dataset{grp, "str", expected.back()}.write(expected, 0, true);
      const auto buff = Dataset{grp, "str"}.read_all<BuffT>();

      REQUIRE(buff.size() == 2);

      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }
    }

    SECTION("iterator") {
      const std::set<std::string> expected{{"a", "b", "c"}};
      Dataset{grp, "str", "a"}.write(expected.begin(), expected.end(), 0, true);
      const auto buff = Dataset{grp, "str"}.read_all<std::vector<std::string>>();

      REQUIRE(buff.size() == 3);

      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(expected.count(buff[i]) == 1);
      }
    }

    SECTION("atomic") {
      const std::string buff = "test";
      Dataset{grp, "str", buff}.write(buff, 3, true);

      CHECK(Dataset{grp, "str"}.read<std::string>(0).empty());
      CHECK(Dataset{grp, "str"}.read<std::string>(3) == buff);
    }
  }

  SECTION("numeric") {
    using T = double;
    using BuffT = std::vector<T>;
    const BuffT expected{0.1, 0.2, 0.3};

    SECTION("vector<N>") {
      Dataset{grp, "num", T{}}.write(expected, 0, true);
      const auto buff = Dataset{grp, "num"}.read_all<BuffT>();

      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }
    }

    SECTION("variant buff") {
      const internal::VariantBuffer vexpected{expected};
      Dataset{grp, "num", T{}}.write(vexpected, 0, true);

      const auto vbuff = Dataset{grp, "num"}.read_all();
      const auto& buff = vbuff.get<T>();

      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }
    }

    SECTION("atomic") {
      Dataset{grp, "num", T{}}.write(7.0, 5, true);
      REQUIRE(Dataset{grp, "num"}.size() == 6);

      CHECK(Dataset{grp, "num"}.read<T>(0) == 0.0);
      CHECK(Dataset{grp, "num"}.read<T>(5) == 7.0);

      const auto vbuff = Dataset{grp, "num"}.read(5);
      CHECK(std::get<T>(vbuff) == 7.0);
    }
  }

  SECTION("out of bound access") {
    CHECK_THROWS_WITH((Dataset{grp, "num", int{}}.write(1, 100U, false)),
                      Catch::Matchers::ContainsSubstring("attempt to access") &&
                          Catch::Matchers::ContainsSubstring("which is empty"));

    Dataset{grp, "num"}.resize(10);
    CHECK_THROWS_WITH((Dataset{grp, "num"}.write(1, 100U, false)),
                      Catch::Matchers::ContainsSubstring("attempt to access") &&
                          Catch::Matchers::ContainsSubstring("past the end"));

    CHECK_THROWS_WITH((Dataset{grp, "num"}.write(std::vector<int>{1, 2, 3}, 100U, false)),
                      Catch::Matchers::ContainsSubstring("attempt to access") &&
                          Catch::Matchers::ContainsSubstring("past the end"));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Dataset - accessors", "[dataset][short]") {
  const auto path = datadir / "cooler_test_file.cool";
  RootGroup grp{HighFive::File(path.string()).getGroup("/")};

  Dataset dset{grp, "chroms/name"};

  CHECK(dset.size() == 20);

  CHECK(dset.get().getFile().getName() == path.string());

  CHECK(dset.file_name() == path.string());
  CHECK(dset.uri() == fmt::format(FMT_STRING("{}::/chroms/name"), path.string()));
  CHECK(dset.hdf5_path() == "/chroms/name");
}

TEST_CASE("Dataset - large read/write", "[dataset][long]") {
  const auto path = testdir() / "test_dataset_large_rw.h5";

  constexpr std::uint64_t seed{4195331987557451569};
  constexpr std::size_t N = 5'000'000;
  {
    RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};
    Dataset dset{grp, "int", std::uint8_t{}};

    std::vector<std::uint8_t> buff(1'000'000);
    std::mt19937_64 rand_eng{seed};
    while (dset.size() != N) {
      std::generate(buff.begin(), buff.end(),
                    [&]() { return static_cast<std::uint8_t>(rand_eng()); });
      dset.append(buff);
    }
    CHECK(dset.size() == N);
  }

  RootGroup grp{HighFive::File(path.string()).getGroup("/")};
  Dataset dset{grp, "int"};
  REQUIRE(dset.size() == N);

  std::mt19937_64 rand_eng{seed};
  auto it = dset.begin<std::uint8_t>();
  for (std::size_t i = 0; i < N; ++i) {
    CHECK(*it++ == static_cast<std::uint8_t>(rand_eng()));
  }
  CHECK(it == dset.end<std::uint8_t>());
}

}  // namespace coolerpp::test::dataset
