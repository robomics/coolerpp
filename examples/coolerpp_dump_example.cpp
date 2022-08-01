// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/std.h>

#include <chrono>
#include <cstdint>
#include <string>

#include "coolerpp/coolerpp.hpp"

using namespace coolerpp;

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
    auto cooler = File::open_read_only(path_to_cooler);
    const auto bin_size = cooler.bin_size();

    using ContactT = std::uint32_t;
    std::for_each(
        cooler.begin<ContactT>(), cooler.end<ContactT>(), [&](const Pixel<ContactT>& pixel) {
          fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:d}\t{}\n"),
                     pixel.coords.chrom1->name, pixel.coords.bin1_start,
                     std::min(pixel.coords.bin1_start + bin_size, pixel.coords.chrom1->size),
                     pixel.coords.chrom2->name, pixel.coords.bin2_start,
                     std::min(pixel.coords.bin2_start + bin_size, pixel.coords.chrom2->size),
                     pixel.count);
        });
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
