// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cassert>
#include <string_view>

#include "coolerpp/bin_table.hpp"
#include "coolerpp/chromosome.hpp"

namespace coolerpp {

constexpr bool PixelCoordinates::operator==(const PixelCoordinates &other) const noexcept {
  return this->_bin1_id == other._bin1_id && this->_bin2_id == other._bin2_id;
}

constexpr bool PixelCoordinates::operator!=(const PixelCoordinates &other) const noexcept {
  return !(*this == other);
}

constexpr bool PixelCoordinates::operator<(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id < other._bin2_id;
  }
  return this->_bin1_id < other._bin1_id;
}

constexpr bool PixelCoordinates::operator<=(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id <= other._bin2_id;
  }
  return this->_bin1_id <= other._bin1_id;
}

constexpr bool PixelCoordinates::operator>(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id > other._bin2_id;
  }
  return this->_bin1_id > other._bin1_id;
}

constexpr bool PixelCoordinates::operator>=(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id >= other._bin2_id;
  }
  return this->_bin1_id >= other._bin1_id;
}

template <typename N>
inline Pixel<N>::operator bool() const noexcept {
  return !!this->coords;
}
template <typename N>
constexpr bool Pixel<N>::operator==(const Pixel<N> &other) const noexcept {
  return this->coords == other.coords && this->count == other.count;
}
template <typename N>
constexpr bool Pixel<N>::operator!=(const Pixel<N> &other) const noexcept {
  return !(*this == other);
}
template <typename N>
constexpr bool Pixel<N>::operator<(const Pixel<N> &other) const noexcept {
  return this->coords < other.coords;
}
template <typename N>
constexpr bool Pixel<N>::operator<=(const Pixel<N> &other) const noexcept {
  return this->coords <= other.coords;
}
template <typename N>
constexpr bool Pixel<N>::operator>(const Pixel<N> &other) const noexcept {
  return this->coords > other.coords;
}
template <typename N>
constexpr bool Pixel<N>::operator>=(const Pixel<N> &other) const noexcept {
  return this->coords >= other.coords;
}

}  // namespace coolerpp

constexpr auto fmt::formatter<coolerpp::PixelCoordinates>::parse(format_parse_context &ctx)
    -> decltype(ctx.begin()) {
  const auto *it = ctx.begin();
  const auto *end = ctx.end();
  const auto fmt_string =
      std::string_view{&(*ctx.begin()), static_cast<std::size_t>(ctx.end() - ctx.begin())};

  if (it != end) {
    if (fmt_string.find("bedpe") != std::string_view::npos) {
      this->presentation = Presentation::bedpe;
      it += std::string_view{"bedpe"}.size();  // NOLINT
    } else if (fmt_string.find("raw") != std::string_view::npos) {
      this->presentation = Presentation::raw;
      it += std::string_view{"raw"}.size();  // NOLINT
    }
  }

  if (it != end && *it != '}') {
    throw fmt::format_error("invalid format");
  }

  return it;
}

template <typename FormatContext>
inline auto fmt::formatter<coolerpp::PixelCoordinates>::format(const coolerpp::PixelCoordinates &c,
                                                               FormatContext &ctx) const
    -> decltype(ctx.out()) {
  if (this->presentation == Presentation::raw) {
    // clang-format off
    return fmt::format_to(ctx.out(),
                          FMT_STRING("{}\t{}"),
                          c.bin1_id(), c.bin2_id());
    // clang-format on
  }

  assert(this->presentation == Presentation::bedpe);
  const auto bin1 = c.bin1();
  const auto bin2 = c.bin2();
  // clang-format off
  return fmt::format_to(ctx.out(),
                        FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}"),
                        bin1.chrom.name,
                        bin1.start,
                        (std::min)(bin1.start + c.bin_size(), bin1.chrom.size),
                        bin2.chrom.name,
                        bin2.start,
                        (std::min)(bin2.start + c.bin_size(), bin2.chrom.size));
  // clang-format on
}

template <typename N>
constexpr auto fmt::formatter<coolerpp::Pixel<N>>::parse(format_parse_context &ctx)
    -> decltype(ctx.begin()) {
  this->coord_formatter.parse(ctx);
  return ctx.end();
}

template <typename N>
constexpr auto fmt::formatter<coolerpp::Pixel<N>>::presentation() const noexcept -> Presentation {
  return this->coord_formatter.presentation;
}

template <typename N>
template <typename FormatContext>
inline auto fmt::formatter<coolerpp::Pixel<N>>::format(const coolerpp::Pixel<N> &p,
                                                       FormatContext &ctx) const
    -> decltype(ctx.out()) {
  if (this->presentation() == Presentation::raw) {
    return fmt::format_to(ctx.out(), FMT_STRING("{:raw}\t{}"), p.coords, p.count);
  }

  assert(this->presentation() == Presentation::bedpe);
  return fmt::format_to(ctx.out(), FMT_STRING("{:bedpe}\t{}"), p.coords, p.count);
}
