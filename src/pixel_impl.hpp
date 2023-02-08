// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cassert>
#include <string_view>

#include "coolerpp/bin_table.hpp"
#include "coolerpp/chromosome.hpp"

namespace coolerpp {

constexpr PixelCoordinates::operator bool() const noexcept {
  return this->chrom1_id != (std::numeric_limits<std::uint32_t>::max)() &&
         this->chrom2_id != (std::numeric_limits<std::uint32_t>::max)();
}

constexpr bool PixelCoordinates::operator==(const PixelCoordinates &other) const noexcept {
  return this->chrom1_id == other.chrom1_id && this->chrom2_id == other.chrom2_id &&
         this->bin1_start == other.bin1_start && this->bin2_start == other.bin2_start;
}

constexpr bool PixelCoordinates::operator!=(const PixelCoordinates &other) const noexcept {
  return !(*this == other);
}

constexpr bool PixelCoordinates::operator<(const PixelCoordinates &other) const noexcept {
  if (this->chrom1_id != other.chrom1_id) {
    return this->chrom1_id < other.chrom1_id;
  }
  if (this->bin1_start != other.bin1_start) {
    return this->bin1_start < other.bin1_start;
  }

  if (this->chrom2_id != other.chrom2_id) {
    return this->chrom2_id < other.chrom2_id;
  }
  return this->bin2_start < other.bin2_start;
}

constexpr bool PixelCoordinates::operator<=(const PixelCoordinates &other) const noexcept {
  if (this->chrom1_id != other.chrom1_id) {
    return this->chrom1_id <= other.chrom1_id;
  }
  if (this->bin1_start != other.bin1_start) {
    return this->bin1_start <= other.bin1_start;
  }

  if (this->chrom2_id != other.chrom2_id) {
    return this->chrom2_id <= other.chrom2_id;
  }
  return this->bin2_start <= other.bin2_start;
}

constexpr bool PixelCoordinates::operator>(const PixelCoordinates &other) const noexcept {
  if (this->chrom1_id != other.chrom1_id) {
    return this->chrom1_id > other.chrom1_id;
  }
  if (this->bin1_start != other.bin1_start) {
    return this->bin1_start > other.bin1_start;
  }

  if (this->chrom2_id != other.chrom2_id) {
    return this->chrom2_id > other.chrom2_id;
  }
  return this->bin2_start > other.bin2_start;
}

constexpr bool PixelCoordinates::operator>=(const PixelCoordinates &other) const noexcept {
  if (this->chrom1_id != other.chrom1_id) {
    return this->chrom1_id >= other.chrom1_id;
  }
  if (this->bin1_start != other.bin1_start) {
    return this->bin1_start >= other.bin1_start;
  }

  if (this->chrom2_id != other.chrom2_id) {
    return this->chrom2_id >= other.chrom2_id;
  }
  return this->bin2_start >= other.bin2_start;
}

template <typename N>
constexpr Pixel<N>::operator bool() const noexcept {
  return !!this->coords;
}
template <typename N>
constexpr bool Pixel<N>::operator==(const Pixel<N> &other) const noexcept {
  return this->coords == other.coords;
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

  // Check if reached the end of the range:
  if (it != end && *it != '}') {
    throw fmt::format_error("invalid format");
  }

  // Return an iterator past the end of the parsed range:
  return it;
}

// Formats the point p using the parsed format specification (presentation)
// stored in this formatter.
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
  // clang-format off
  return fmt::format_to(ctx.out(),
                        FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}"),
                        c.chrom1().name,
                        c.bin1_start,
                        (std::min)(c.bin1_start + c.bin_size(), c.chrom1().size),
                        c.chrom2().name,
                        c.bin2_start,
                        (std::min)(c.bin2_start + c.bin_size(), c.chrom2().size));
  // clang-format on
}
