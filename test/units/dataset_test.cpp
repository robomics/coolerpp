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
TEST_CASE("Dataset: read", "[dataset][short]") {
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
      std::uint64_t buff{};
      dset.read(buff, 2);
      CHECK(buff == 159'599'783);

      CHECK(dset.read_last<std::int32_t>() == 166'650'296);
      CHECK(std::get<std::int32_t>(dset.read_last()) == 166'650'296);
    }

    SECTION("enum") { CHECK(Dataset{grp, "bins/chrom"}.read<std::uint32_t>(0) == 0); }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Dataset: write", "[dataset][short]") {
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

      for (const auto& i : buff) {
        CHECK(expected.count(i) == 1);
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
TEST_CASE("Dataset: accessors", "[dataset][short]") {
  const auto path = datadir / "cooler_test_file.cool";
  RootGroup grp{HighFive::File(path.string()).getGroup("/")};

  Dataset dset{grp, "chroms/name"};

  CHECK(dset.size() == 20);

  CHECK(dset.get().getFile().getName() == path.string());

  CHECK(dset.file_name() == path.string());
  CHECK(dset.uri() == fmt::format(FMT_STRING("{}::/chroms/name"), path.string()));
  CHECK(dset.hdf5_path() == "/chroms/name");
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Dataset: linear iteration", "[dataset][short]") {
  const auto path1 = datadir / "cooler_test_file.cool";

  RootGroup grp1{HighFive::File(path1.string()).getGroup("/")};
  Dataset dset1(grp1, "/pixels/count");

  std::vector<std::uint32_t> pixel_buff;
  dset1.read_all(pixel_buff);
  REQUIRE(pixel_buff.size() == 107'041);

  SECTION("forward") {
    auto it = dset1.begin<std::uint32_t>();
    auto last_pixel = dset1.end<std::uint32_t>();
    REQUIRE(std::distance(it, last_pixel) == 107'041);

    for (const auto& expected : pixel_buff) {
      REQUIRE(it != last_pixel);
      CHECK(*it++ == expected);
    }
    CHECK(it == last_pixel);
  }

  SECTION("backward") {
    auto it = dset1.end<std::uint32_t>();
    auto first_pixel = dset1.begin<std::uint32_t>();

    REQUIRE(std::distance(first_pixel, it) == 107'041);

    for (std::size_t i = pixel_buff.size(); i != 0; --i) {
      REQUIRE(it != first_pixel);
      CHECK(*(--it) == pixel_buff[i - 1]);
    }

    CHECK(it == first_pixel);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Dataset: random iteration", "[dataset][medium]") {
  const auto path = testdir() / "dataset_iterator_random.h5";

  RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};
  Dataset dset(grp, "int", std::uint64_t{});

  std::random_device rd;
  std::mt19937_64 rand_eng{rd()};

  constexpr std::size_t N = 5'000'000;
  std::vector<std::uint64_t> buff(N);
  std::generate(buff.begin(), buff.end(), [&]() { return rand_eng(); });
  dset.append(buff);
  REQUIRE(dset.size() == N);

  SECTION("operator -/+") {
    auto first = dset.begin<std::uint64_t>();
    auto last = dset.end<std::uint64_t>();
    for (std::size_t i = 0; i < 100; ++i) {
      const auto j = std::uniform_int_distribution<std::uint64_t>{0, N - 1}(rand_eng);

      CHECK(*(first + j) == buff[j]);
      CHECK(*(last - j) == buff[N - j]);
    }
  }

  SECTION("subsequent calls to operator+=") {
    for (std::size_t i = 0; i < 10; ++i) {
      auto first = dset.begin<std::uint64_t>();
      auto last = dset.end<std::uint64_t>();
      std::size_t j = 0;

      while (first < last) {
        CHECK(*first == buff[j]);

        const auto step = std::uniform_int_distribution<std::size_t>{
            0, std::min(std::size_t(500), buff.size() - j)}(rand_eng);
        first += step;
        j += step;
      }
    }
  }

  SECTION("subsequent calls to operator-=") {
    for (std::size_t i = 0; i < 10; ++i) {
      auto first = dset.end<std::uint64_t>() - 1;
      auto last = dset.begin<std::uint64_t>();
      std::size_t j = buff.size() - 1;

      while (first > last) {
        CHECK(*first == buff[j]);

        const auto step =
            std::uniform_int_distribution<std::size_t>{0, std::min(std::size_t(500), j)}(rand_eng);

        first -= step;
        j -= step;
      }
    }
  }
}

TEST_CASE("Dataset: large read/write", "[dataset][long]") {
  const auto path = testdir() / "test_dataset_large_rw.h5";

  constexpr std::uint64_t seed{4195331987557451569};
  constexpr std::size_t N = 5'000'000;
  {
    RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};
    Dataset dset{grp, "int", std::uint8_t{}};

    std::vector<std::uint8_t> buff(1'000'000);
    std::mt19937_64 rand_eng{seed};  // NOLINT(cert-msc32-c,cert-msc51-cpp)
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

  std::mt19937_64 rand_eng{seed};  // NOLINT(cert-msc32-c,cert-msc51-cpp)
  auto it = dset.begin<std::uint8_t>();
  for (std::size_t i = 0; i < N; ++i) {
    CHECK(*it++ == static_cast<std::uint8_t>(rand_eng()));
  }
  CHECK(it == dset.end<std::uint8_t>());
}

}  // namespace coolerpp::test::dataset
