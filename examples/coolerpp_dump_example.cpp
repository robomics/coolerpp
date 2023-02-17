// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/os.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <string>

#include "coolerpp/coolerpp.hpp"

using namespace coolerpp;

struct Coords {
  std::string chrom{};
  std::uint32_t start{};
  std::uint32_t end{};
};

template <typename It>
static std::size_t print_pixels(It first_pixel, It last_pixel) {
  std::size_t nnz = 0;
  std::for_each(first_pixel, last_pixel, [&](const auto& pixel) {
    ++nnz;
    fmt::print(FMT_COMPILE("{:bedpe}\t{}\n"), pixel.coords, pixel.count);
  });

  return nnz;
}

static std::size_t dump(const File& cooler, const std::string& coord1, const std::string& coord2) {
  if (cooler.has_integral_pixels()) {
    auto selector = cooler.fetch<std::int64_t>(coord1, coord2);
    return print_pixels(selector.begin(), selector.end());
  }

  auto selector = cooler.fetch<double>(coord1, coord2);
  return print_pixels(selector.begin(), selector.end());
}

static std::size_t dump(const File& cooler, const std::string& coords) {
  return dump(cooler, coords, coords);
}

static std::size_t dump(const File& cooler) {
  if (cooler.has_integral_pixels()) {
    return print_pixels(cooler.begin<std::int64_t>(), cooler.end<std::int64_t>());
  }
  return print_pixels(cooler.begin<double>(), cooler.end<double>());
}

int main(int argc, char** argv) {
  if (argc < 2) {
    fmt::print(stderr,
               FMT_STRING("Usage:   {0} my_cooler.cool [region1] [region2]\n"
                          "Example: {0} my_cooler.cool\n"
                          "Example: {0} my_cooler.mcool::/resolutions/10000\n"
                          "Example: {0} my_cooler.cool chr1\n"
                          "Example: {0} my_cooler.cool chr1 chr2\n"
                          "Example: {0} my_cooler.cool chr1:50000-100000\n"
                          "Example: {0} my_cooler.cool chr1:50000-100000 chr2\n"),
               argv[0]);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    return 1;
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const std::string path_to_cooler = argv[1];

  try {
    const auto t0 = std::chrono::steady_clock::now();
    const auto cooler = File::open_read_only(path_to_cooler);
    std::size_t nnz = 0;

    if (argc == 2) {
      nnz = dump(cooler);
    } else if (argc == 3) {
      nnz = dump(cooler, argv[2]);
    } else {
      assert(argc == 4);
      nnz = dump(cooler, argv[2], argv[3]);
    }

    const auto t1 = std::chrono::steady_clock::now();

    const auto elapsed_time_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    fmt::print(stderr, FMT_STRING("Dumped {} pixels in {}s!\n"), nnz,
               static_cast<double>(elapsed_time_ms) / 1000.0);
  } catch (const std::exception& e) {
    fmt::print(
        stderr,
        FMT_STRING("The following error occurred while running coolerpp_dump on file {}: {}\n"),
        path_to_cooler, e.what());
    return 1;
  }
}
