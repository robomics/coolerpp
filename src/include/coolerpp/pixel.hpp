// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstdint>
#include <string_view>
#include <type_traits>

namespace coolerpp {

struct Chromosome;
class BinTableLazy;

class PixelCoordinates {
  const BinTableLazy *_bins{};

 public:
  const Chromosome *chrom1{};
  const Chromosome *chrom2{};

  std::uint32_t bin1_start{};
  std::uint32_t bin2_start{};

  PixelCoordinates() = delete;
  PixelCoordinates(const BinTableLazy &bins, const Chromosome &chrom1_, const Chromosome &chrom2_,
                   std::uint32_t bin1_start_, std::uint32_t bin2_start_);
  PixelCoordinates(const BinTableLazy &bins, std::string_view chrom1_name,
                   std::string_view chrom2_name, std::uint32_t bin1_start_,
                   std::uint32_t bin2_start_);
  PixelCoordinates(const BinTableLazy &bins, std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                   std::uint32_t bin1_start_, std::uint32_t bin2_start_);

  PixelCoordinates(const BinTableLazy &bins, const Chromosome &chrom, std::uint32_t bin1_start_,
                   std::uint32_t bin2_start_);
  PixelCoordinates(const BinTableLazy &bins, std::string_view chrom_name, std::uint32_t bin1_start_,
                   std::uint32_t bin2_start_);
  PixelCoordinates(const BinTableLazy &bins, std::uint32_t chrom_id, std::uint32_t bin1_start_,
                   std::uint32_t bin2_start_);

  PixelCoordinates(const BinTableLazy &bins, std::uint64_t bin1_id, std::uint64_t bin2_id);

  [[nodiscard]] bool operator==(const PixelCoordinates &other) const;
  [[nodiscard]] bool operator<(const PixelCoordinates &other) const;

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

  [[nodiscard]] inline bool operator==(const Pixel &other) const {
    return this->coords == other.coords && this->count == other.count;
  }
  [[nodiscard]] inline bool operator<(const Pixel &other) const {
    return this->coords < other.coords;
  }
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
