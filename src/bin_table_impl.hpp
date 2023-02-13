// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

namespace coolerpp {  // NOLINT

constexpr Bin::Bin(const Chromosome &chrom_, std::uint32_t start_, std::uint32_t end_) noexcept
    : chrom(chrom_), start(start_), end(end_) {
  assert(start <= end);
}

template <class ChromIt>
inline BinTableLazy::BinTableLazy(ChromIt first_chrom, ChromIt last_chrom, std::uint32_t bin_size)
    : BinTableLazy(ChromosomeSet(first_chrom, last_chrom), bin_size) {}

template <class ChromNameIt, class ChromSizeIt>
inline BinTableLazy::BinTableLazy(ChromNameIt first_chrom_name, ChromNameIt last_chrom_name,
                                  ChromSizeIt first_chrom_size, std::uint32_t bin_size)
    : BinTableLazy(ChromosomeSet(first_chrom_name, last_chrom_name, first_chrom_size), bin_size) {}

constexpr std::uint32_t BinTableLazy::bin_size() const noexcept { return this->_bin_size; }

constexpr const ChromosomeSet &BinTableLazy::chromosomes() const noexcept { return this->_chroms; }

constexpr const std::vector<std::uint64_t> &BinTableLazy::num_bin_prefix_sum() const noexcept {
  return this->_num_bins_prefix_sum;
}

constexpr auto BinTableLazy::begin() -> iterator { return iterator(*this); }
constexpr auto BinTableLazy::end() -> iterator { return iterator::make_end_iterator(*this); }
constexpr auto BinTableLazy::begin() const -> const_iterator { return const_iterator(*this); }
constexpr auto BinTableLazy::end() const -> const_iterator {
  return const_iterator::make_end_iterator(*this);
}
constexpr auto BinTableLazy::cbegin() const -> const_iterator { return this->begin(); }
constexpr auto BinTableLazy::cend() const -> const_iterator { return this->end(); }

constexpr BinTableLazy::iterator::iterator(const BinTableLazy &bin_table) noexcept
    : _bin_table{&bin_table} {}

constexpr bool BinTableLazy::iterator::operator==(const iterator &other) const noexcept {
  // clang-format off
  return this->_bin_table == other._bin_table &&
         this->_chrom_id == other._chrom_id &&
         this->_idx == other._idx;
  // clang-format on
}

constexpr bool BinTableLazy::iterator::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

constexpr auto BinTableLazy::iterator::make_end_iterator(const BinTableLazy &table) noexcept
    -> iterator {
  iterator it(table);

  it._chrom_id = std::numeric_limits<std::uint32_t>::max();
  it._idx = npos;
  return it;
}

constexpr std::uint32_t BinTableLazy::iterator::bin_size() const noexcept {
  assert(this->_bin_table);
  return this->_bin_table->bin_size();
}

}  // namespace coolerpp

constexpr auto fmt::formatter<coolerpp::Bin>::parse(format_parse_context &ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.end();
}

template <class FormatContext>
inline auto fmt::formatter<coolerpp::Bin>::format(const coolerpp::Bin &b, FormatContext &ctx) const
    -> decltype(ctx.out()) {
  return fmt::format_to(ctx.out(), FMT_STRING("{}:{}-{}"), b.chrom.name, b.start, b.end);
}
