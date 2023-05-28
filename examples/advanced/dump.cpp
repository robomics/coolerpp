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

static void print(const Chromosome& chrom) {
  fmt::print(FMT_COMPILE("{:s}\t{:d}\n"), chrom.name, chrom.size);
}

static void print(const Bin& bin) {
  fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\n"), bin.chrom.name, bin.start, bin.end);
}

template <bool join>
static void print(const Pixel<std::int64_t>& pixel) {
  if constexpr (join) {
    fmt::print(FMT_COMPILE("{:bedpe}\t{:d}\n"), pixel.coords, pixel.count);
  } else {
    fmt::print(FMT_COMPILE("{:raw}\t{:d}\n"), pixel.coords, pixel.count);
  }
}

template <bool join>
static void print(const Pixel<double>& pixel) {
  if constexpr (join) {
    fmt::print(FMT_COMPILE("{:bedpe}\t{:g}\n"), pixel.coords, pixel.count);
  } else {
    fmt::print(FMT_COMPILE("{:raw}\t{:g}\n"), pixel.coords, pixel.count);
  }
}

static void dump_chroms(const File& clr, std::string_view range) {
  if (range == "all") {
    for (const auto& chrom : clr.chromosomes()) {
      print(chrom);
    }
    return;
  }

  const auto coords = PixelSelector<int>::parse_query(clr.bins_ptr(), range);
  auto it = clr.chromosomes().find(coords.chrom1());
  if (it != clr.chromosomes().end()) {
    print(*it);
  }
}

static void dump_bins(const File& clr, std::string_view range) {
  if (range == "all") {
    for (const auto& bin : clr.bins()) {
      print(bin);
    }
    return;
  }

  const auto coords = PixelSelector<int>::parse_query(clr.bins_ptr(), range);
  const auto bins = clr.bins().at(coords.chrom1());

  auto first_bin = std::lower_bound(bins.begin(), bins.end(), coords.bin1());
  auto last_bin = std::upper_bound(first_bin, bins.end(), coords.bin2());

  std::for_each(first_bin, last_bin, [](const Bin& bin) { print(bin); });
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
