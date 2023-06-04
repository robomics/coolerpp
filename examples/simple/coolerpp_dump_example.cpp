// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <exception>
#include <string_view>

#include "coolerpp/coolerpp.hpp"

using namespace coolerpp;

template <typename It>
static void print_pixels(It first_pixel, It last_pixel) {
  std::for_each(first_pixel, last_pixel,
                [](const auto& pixel) { fmt::print(FMT_STRING("{:bg2}\n"), pixel); });
}

template <typename N>
static void dump(const File& cooler, std::string_view range1, std::string_view range2) {
  auto selector = cooler.fetch<N>(range1, range2);
  print_pixels(selector.begin(), selector.end());
}

template <typename N>
static void dump(const File& cooler) {
  print_pixels(cooler.begin<N>(), cooler.end<N>());
}

static void print_usage(std::string_view arg0) {
  fmt::print(stderr,
             FMT_STRING("Usage:   {0} my_cooler.cool [region1] [region2]\n"
                        "Example: {0} my_cooler.cool\n"
                        "Example: {0} my_cooler.mcool::/resolutions/10000\n"
                        "Example: {0} my_cooler.cool chr1\n"
                        "Example: {0} my_cooler.cool chr1 chr2\n"
                        "Example: {0} my_cooler.cool chr1:50000-100000\n"
                        "Example: {0} my_cooler.cool chr1:50000-100000 chr2\n"),
             arg0);
}

int main(int argc, char** argv) {
  if (argc < 2) {
    print_usage(argv[0]);  // NOLINT
    return 1;
  }

  const std::string_view path_to_cooler = argv[1];              // NOLINT
  const std::string_view range1 = argc < 3 ? "" : argv[2];      // NOLINT;
  const std::string_view range2 = argc < 4 ? range1 : argv[3];  // NOLINT;

  if (path_to_cooler == "--help" || path_to_cooler == "-h") {
    print_usage(argv[0]);  // NOLINT
    return 0;
  }

  try {
    const auto cooler = File::open_read_only(path_to_cooler);
    const auto has_int_pixels = cooler.has_integral_pixels();

    if (range1.empty() && range2.empty()) {
      has_int_pixels ? dump<std::int64_t>(cooler) : dump<double>(cooler);
    } else {
      has_int_pixels ? dump<std::int64_t>(cooler, range1, range2)
                     : dump<double>(cooler, range1, range2);
    }
  } catch (const std::exception& e) {
    fmt::print(
        stderr,
        FMT_STRING("The following error occurred while running coolerpp_dump on URI \"{}\": {}\n"),
        path_to_cooler, e.what());
    return 1;
  }
}
