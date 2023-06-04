// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <string_view>

#include "coolerpp/coolerpp.hpp"
#include "coolerpp_tools/config.hpp"
#include "coolerpp_tools/tools.hpp"

namespace coolerpp::tools {

template <bool join>
static void print(const Pixel<std::int64_t>& pixel) {
  if constexpr (join) {
    fmt::print(FMT_COMPILE("{:bg2}\n"), pixel);
  } else {
    fmt::print(FMT_COMPILE("{:raw}\n"), pixel);
  }
}

template <bool join>
static void print(const Pixel<double>& pixel) {
  if constexpr (join) {
    fmt::print(FMT_COMPILE("{:bg2}\t{:g}\n"), pixel.coords, pixel.count);
  } else {
    fmt::print(FMT_COMPILE("{:raw}\t{:g}\n"), pixel.coords, pixel.count);
  }
}

static void dump_chroms(const File& clr, std::string_view range) {
  if (range == "all") {
    for (const Chromosome& chrom : clr.chromosomes()) {
      fmt::print(FMT_COMPILE("{:s}\t{:d}\n"), chrom.name(), chrom.size());
    }
    return;
  }

  const auto coords = GenomicInterval::parse_ucsc(clr.chromosomes(), std::string{range});
  auto it = clr.chromosomes().find(coords.chrom());
  if (it != clr.chromosomes().end()) {
    fmt::print(FMT_COMPILE("{:s}\t{:d}\n"), it->name(), it->size());
  }
}

static void dump_bins(const File& clr, std::string_view range) {
  if (range == "all") {
    for (const auto& bin : clr.bins()) {
      fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\n"), bin.chrom().name(), bin.start(), bin.end());
    }
    return;
  }

  const auto coords = GenomicInterval::parse_ucsc(clr.chromosomes(), std::string{range});
  auto [first_bin, last_bin] = clr.bins().find_overlap(coords);
  std::for_each(first_bin, last_bin, [](const Bin& bin) {
    fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\n"), bin.chrom().name(), bin.start(), bin.end());
  });
}

template <typename N, bool join>
static void print_pixels(typename PixelSelector<N>::iterator first_pixel,
                         typename PixelSelector<N>::iterator last_pixel,
                         const std::shared_ptr<const Weights>& weights) {
  if (weights != nullptr) {
    const auto sel = Balancer<N>(first_pixel, last_pixel, weights);
    std::for_each(sel.begin(), sel.end(), [&](const auto& pixel) { print<join>(pixel); });
    return;
  }
  std::for_each(first_pixel, last_pixel, [&](const auto& pixel) { print<join>(pixel); });
}

template <bool join>
static void dump_pixels(const File& clr, std::string_view range1, std::string_view range2,
                        std::string_view balanced) {
  const auto has_int_pixels = clr.has_integral_pixels();
  const auto weights =
      balanced.empty() ? std::shared_ptr<const Weights>(nullptr) : clr.read_weights(balanced);

  if (range1 == "all") {
    assert(range2 == "all");
    return has_int_pixels
               ? print_pixels<std::int64_t, join>(clr.begin<std::int64_t>(),
                                                  clr.end<std::int64_t>(), weights)
               : print_pixels<double, join>(clr.begin<double>(), clr.end<double>(), weights);
  }

  if (has_int_pixels) {
    auto sel = clr.fetch<std::int64_t>(range1, range2);
    return print_pixels<std::int64_t, join>(sel.begin(), sel.end(), weights);
  }

  auto sel = clr.fetch<double>(range1, range2);
  return print_pixels<double, join>(sel.begin(), sel.end(), weights);
}

void dump_subcmd(const DumpConfig& c) {
  const auto clr = File::open_read_only(c.uri);

  if (c.table == "chroms") {
    return dump_chroms(clr, c.range1);
  }
  if (c.table == "bins") {
    return dump_bins(clr, c.range1);
  }

  assert(c.table == "pixels");
  if (c.join) {
    dump_pixels<true>(clr, c.range1, c.range2, c.balanced);
  } else {
    dump_pixels<false>(clr, c.range1, c.range2, c.balanced);
  }
}

}  // namespace coolerpp::tools
