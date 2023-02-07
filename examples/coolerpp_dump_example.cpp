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

template <class N>
void print_pixel(const Pixel<N>& pixel) {
  fmt::print(FMT_COMPILE("{:bedpe}\t{}\n"), pixel.coords, pixel.count);
}

int main(int argc, char** argv) {
  if (argc != 2) {
    fmt::print(stderr,
               FMT_STRING("Usage:   {0} my_cooler.cool\n"
                          "Example: {0} my_cooler.cool\n"
                          "Example: {0} my_cooler.mcool::/resolutions/10000\n"),
               argv[0]);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    return 1;
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const std::string path_to_cooler = argv[1];

  try {
    const auto t0 = std::chrono::steady_clock::now();
    const auto cooler = File::open_read_only(path_to_cooler);

    if (cooler.has_integral_pixels()) {
      using T = std::int64_t;
      std::for_each(cooler.begin<T>(), cooler.end<T>(), print_pixel<T>);

    } else {
      using T = double;
      std::for_each(cooler.begin<T>(), cooler.end<T>(), print_pixel<T>);
    }

    const auto t1 = std::chrono::steady_clock::now();

    const auto elapsed_time_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    fmt::print(stderr, FMT_STRING("Dumped {} pixels in {}s!\n"), *cooler.attributes().nnz,
               static_cast<double>(elapsed_time_ms) / 1000.0);
  } catch (const std::exception& e) {
    fmt::print(
        stderr,
        FMT_STRING("The following error occurred while running coolerpp_dump on file {}: \"{}\"\n"),
        path_to_cooler, e.what());
    return 1;
  }
}
