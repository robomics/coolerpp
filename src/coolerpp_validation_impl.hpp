// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <highfive/H5Exception.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

#include "coolerpp/utils.hpp"

namespace coolerpp {

inline void File::validate_bins() const {
  try {
    assert(this->_attrs.bin_type == "fixed");
    auto nchroms = this->dataset("bins/chrom").size();
    auto nstarts = this->dataset("bins/start").size();
    auto nends = this->dataset("bins/end").size();
    if (nchroms != nstarts || nchroms != nends) {
      throw std::runtime_error(fmt::format(FMT_STRING("Datasets have inconsistent sizes:\n"
                                                      " - \"bins/chrom\": {}\n"
                                                      " - \"bins/start\": {}\n"
                                                      " - \"bins/end\": {}\n"
                                                      "Expected {}"),
                                           nchroms, nstarts, nends, this->bins().size()));
    }

    const auto &nbins = nchroms;
    if (nbins != this->bins().size()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Expected {} bins, found {}"), this->bins().size(), nchroms));
    }

    auto chrom_it = this->dataset("bins/chrom").begin<std::uint32_t>();
    auto start_it = this->dataset("bins/start").begin<std::uint32_t>();
    auto end_it = this->dataset("bins/end").begin<std::uint32_t>();

    auto last_chrom = this->dataset("bins/chrom").end<std::uint32_t>();
    auto last_start = this->dataset("bins/start").end<std::uint32_t>();
    auto last_end = this->dataset("bins/end").end<std::uint32_t>();

    std::size_t i = 0;
    for (const Bin &bin : this->bins()) {
      if (chrom_it == last_chrom || start_it == last_start || end_it == last_end) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Expected {} bins, found {}"), this->bins().size(), i));
      }

      if (this->chromosomes().at(*chrom_it).name != bin.chrom.name || *start_it != bin.start ||
          *end_it != bin.end) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Bin #{}: expected {}:{}-{}, found {}:{}-{}"), i,
                        this->chromosomes().at(*chrom_it).name, *start_it, *end_it, bin.chrom.name,
                        bin.start, bin.end));
      }
      ++chrom_it;
      ++start_it;
      ++end_it;
      ++i;
    }

  } catch (const HighFive::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Bin table at URI {}/{} is invalid or corrupted: {}"), this->uri(),
                    this->group("bins")().getPath(), e.what()));
  }
}

template <typename PixelIt>
inline void File::validate_pixels_before_append(PixelIt first_pixel, PixelIt last_pixel) const {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using T = decltype(std::declval<PixelT>().count);
  try {
    std::for_each(first_pixel, last_pixel, [&](const Pixel<T> &pixel) {
      if (pixel.count == T{0}) {
        throw std::runtime_error("found a pixel of value 0");
      }

      if (!this->chromosomes().contains(pixel.coords.chrom1_id())) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("invalid chromosome id {}"), pixel.coords.chrom1_id()));
      }

      if (pixel.coords.chrom1_id() != pixel.coords.chrom2_id() &&
          !this->chromosomes().contains(pixel.coords.chrom2_id())) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("invalid chromosome id {}"), pixel.coords.chrom2_id()));
      }

      if (const auto bin_id = pixel.coords.bin1_id(); bin_id > this->bins().size()) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("invalid bin id {}: bin maps outside of the bin table"), bin_id));
      }

      if (const auto bin_id = pixel.coords.bin2_id(); bin_id > this->bins().size()) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("invalid bin id {}: bin maps outside of the bin table"), bin_id));
      }

      if (pixel.coords.bin1_id() > pixel.coords.bin2_id()) {
        throw std::runtime_error(fmt::format(FMT_STRING("bin1_id is greater than bin2_id: {} > {}"),
                                             pixel.coords.bin1_id(), pixel.coords.bin2_id()));
      }
    });

    if (!this->dataset("pixels/bin1_id").empty()) {
      const auto last_bin1 = this->dataset("pixels/bin1_id").read_last<std::size_t>();
      const auto last_bin2 = this->dataset("pixels/bin2_id").read_last<std::size_t>();

      const auto new_bin1 = first_pixel->coords.bin1_id();
      const auto new_bin2 = first_pixel->coords.bin2_id();

      if (last_bin1 == new_bin1) {
        if (last_bin2 >= new_bin2) {
          const auto coord1 = this->bins().bin_id_to_coords(new_bin2);
          const auto coord2 = this->bins().bin_id_to_coords(last_bin2);
          throw std::runtime_error(fmt::format(
              FMT_STRING("new pixel {} is located upstream of pixel {}"), coord1, coord2));
        }
      } else if (last_bin1 >= new_bin1) {
        const auto coord1 = this->bins().bin_id_to_coords(new_bin1);
        const auto coord2 = this->bins().bin_id_to_coords(last_bin1);
        throw std::runtime_error(fmt::format(
            FMT_STRING("new pixel {} is located upstream of pixel {}"), coord1, coord2));
      }
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("pixel validation failed: {}"), e.what()));
  }
}

template <typename PixelT>
inline void File::validate_pixel_type() const noexcept {
  static_assert(std::is_arithmetic_v<PixelT>);

  auto assert_holds_alternative = [](const auto &buff, [[maybe_unused]] auto alt) {
    using T [[maybe_unused]] = decltype(alt);
    if (buff.has_value()) {
      assert(std::holds_alternative<T>(*buff));
    }
  };

  if constexpr (std::is_floating_point_v<PixelT>) {
    assert(this->has_float_pixels());
    assert_holds_alternative(this->_attrs.sum, double{});
    assert_holds_alternative(this->_attrs.cis, double{});
  } else if constexpr (std::is_signed_v<PixelT>) {
    assert(this->has_signed_pixels());
    assert_holds_alternative(this->_attrs.sum, std::int64_t{});
    assert_holds_alternative(this->_attrs.cis, std::int64_t{});
  } else {
    assert(this->has_unsigned_pixels());
    assert_holds_alternative(this->_attrs.sum, std::uint64_t{});
    assert_holds_alternative(this->_attrs.cis, std::uint64_t{});
  }
}

}  // namespace coolerpp
