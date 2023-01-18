// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cassert>
#include <string_view>

#include "coolerpp/bin_table.hpp"
#include "coolerpp/chromosome.hpp"

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
  assert(c.chrom1);
  assert(c.chrom2);

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
                        c.chrom1->name,
                        c.bin1_start,
                        std::min(c.bin1_start + c.bin_size(), c.chrom1->size),
                        c.chrom2->name,
                        c.bin2_start,
                        std::min(c.bin2_start + c.bin_size(), c.chrom2->size));
  // clang-format on
}