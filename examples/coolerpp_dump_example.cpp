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

template <typename N>
static std::size_t print_pixels(typename PixelSelector<N>::iterator first_pixel,
                                typename PixelSelector<N>::iterator last_pixel,
                                const std::shared_ptr<Weights>& weights) {
  std::size_t nnz = 0;
  if (weights != nullptr) {
    const auto sel = Balancer<N>(first_pixel, last_pixel, weights);
    std::for_each(sel.begin(), sel.end(), [&](const auto& pixel) {
      ++nnz;
      fmt::print(FMT_COMPILE("{:bg2}\t{}\n"), pixel.coords, pixel.count);
    });
  } else {
    std::for_each(first_pixel, last_pixel, [&](const auto& pixel) {
      ++nnz;
      fmt::print(FMT_COMPILE("{:bg2}\t{}\n"), pixel.coords, pixel.count);
    });
  }

  return nnz;
}

static std::size_t dump(const File& cooler, const std::string& coord1, const std::string& coord2,
                        const std::shared_ptr<Weights>& weights) {
  if (cooler.has_integral_pixels()) {
    auto selector = cooler.fetch<std::int64_t>(coord1, coord2);
    return print_pixels<std::int64_t>(selector.begin(), selector.end(), weights);
  }

  auto selector = cooler.fetch<double>(coord1, coord2);
  return print_pixels<double>(selector.begin(), selector.end(), weights);
}

static std::size_t dump(const File& cooler, const std::shared_ptr<Weights>& weights) {
  if (cooler.has_integral_pixels()) {
    return print_pixels<std::int64_t>(cooler.begin<std::int64_t>(), cooler.end<std::int64_t>(),
                                      weights);
  }
  return print_pixels<double>(cooler.begin<double>(), cooler.end<double>(), weights);
}

int main(int argc, char** argv) {
  if (argc < 3) {
    fmt::print(stderr,
               FMT_STRING("Usage:   {0} my_cooler.cool balancing [region1] [region2]\n"
                          "Example: {0} my_cooler.cool raw\n"
                          "Example: {0} my_cooler.cool weight\n"
                          "Example: {0} my_cooler.mcool::/resolutions/10000 raw\n"
                          "Example: {0} my_cooler.cool raw chr1\n"
                          "Example: {0} my_cooler.cool raw chr1 chr2\n"
                          "Example: {0} my_cooler.cool raw chr1:50000-100000\n"
                          "Example: {0} my_cooler.cool raw chr1:50000-100000 chr2\n"),
               argv[0]);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    return 1;
  }

  const std::string path_to_cooler = argv[1];  // NOLINT
  const std::string balancing = argv[2];       // NOLINT

  try {
    const auto t0 = std::chrono::steady_clock::now();
    const auto clr = File::open_read_only_read_once(path_to_cooler);
    const auto weights =
        balancing == "raw" ? std::shared_ptr<Weights>(nullptr) : clr.read_weights(balancing);
    std::size_t nnz = 0;

    if (argc == 3) {
      nnz = dump(clr, weights);
    } else if (argc == 4) {
      nnz = dump(clr, argv[3], argv[3], weights);  // NOLINT
    } else {
      assert(argc == 5);
      nnz = dump(clr, argv[3], argv[4], weights);  // NOLINT
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
