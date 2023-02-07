// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstdint>
#include <limits>
#include <string_view>
#include <type_traits>

namespace coolerpp {

struct Chromosome;
class BinTableLazy;

class PixelCoordinates {
  const BinTableLazy *_bins{};

 public:
  std::uint32_t chrom1_id{(std::numeric_limits<std::uint32_t>::max)()};  // NOLINT
  std::uint32_t chrom2_id{(std::numeric_limits<std::uint32_t>::max)()};  // NOLINT

  std::uint32_t bin1_start{};  // NOLINT
  std::uint32_t bin2_start{};  // NOLINT

  PixelCoordinates() = delete;
  PixelCoordinates(const BinTableLazy &bins, const Chromosome &chrom1, const Chromosome &chrom2,
                   std::uint32_t bin1_start_, std::uint32_t bin2_start_);
  PixelCoordinates(const BinTableLazy &bins, std::string_view chrom1_name,
                   std::string_view chrom2_name, std::uint32_t bin1_start_,
                   std::uint32_t bin2_start_);
  PixelCoordinates(const BinTableLazy &bins, std::uint32_t chrom1_id_, std::uint32_t chrom2_id_,
                   std::uint32_t bin1_start_, std::uint32_t bin2_start_);

  PixelCoordinates(const BinTableLazy &bins, const Chromosome &chrom, std::uint32_t bin1_start_,
                   std::uint32_t bin2_start_);
  PixelCoordinates(const BinTableLazy &bins, std::string_view chrom_name, std::uint32_t bin1_start_,
                   std::uint32_t bin2_start_);
  PixelCoordinates(const BinTableLazy &bins, std::uint32_t chrom_id_, std::uint32_t bin1_start_,
                   std::uint32_t bin2_start_);

  PixelCoordinates(const BinTableLazy &bins, std::uint64_t bin1_id, std::uint64_t bin2_id);

  [[nodiscard]] constexpr explicit operator bool() const noexcept;
  [[nodiscard]] constexpr bool operator==(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] constexpr bool operator!=(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] constexpr bool operator<(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] constexpr bool operator<=(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] constexpr bool operator>(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] constexpr bool operator>=(const PixelCoordinates &other) const noexcept;

  [[nodiscard]] const Chromosome &chrom1() const;
  [[nodiscard]] const Chromosome &chrom2() const;

  [[nodiscard]] std::uint64_t bin1_id() const;
  [[nodiscard]] std::uint64_t bin2_id() const;

  [[nodiscard]] std::uint32_t bin_size() const noexcept;
};

template <class N>
struct Pixel {
  static_assert(std::is_arithmetic_v<N>);
  using Coordinates = PixelCoordinates;

  Coordinates coords;
  N count;

  [[nodiscard]] constexpr explicit operator bool() const noexcept;
  [[nodiscard]] constexpr bool operator==(const Pixel<N> &other) const noexcept;
  [[nodiscard]] constexpr bool operator!=(const Pixel<N> &other) const noexcept;
  [[nodiscard]] constexpr bool operator<(const Pixel<N> &other) const noexcept;
  [[nodiscard]] constexpr bool operator<=(const Pixel<N> &other) const noexcept;
  [[nodiscard]] constexpr bool operator>(const Pixel<N> &other) const noexcept;
  [[nodiscard]] constexpr bool operator>=(const Pixel<N> &other) const noexcept;
};

}  // namespace coolerpp

template <>
struct fmt::formatter<coolerpp::PixelCoordinates> {
  // Presentation can be any of the following:
  //  - raw
  //  - bedpe

  constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin());
  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  auto format(const coolerpp::PixelCoordinates &c, FormatContext &ctx) const -> decltype(ctx.out());

 private:
  enum Presentation { raw, bedpe };
  Presentation presentation{Presentation::bedpe};
};

#include "../../pixel_impl.hpp"
