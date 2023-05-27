// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "coolerpp/common.hpp"

namespace coolerpp {  // NOLINT

constexpr Bin::Bin(const Chromosome &chrom_, std::uint32_t start_, std::uint32_t end_) noexcept
    : chrom(chrom_), start(start_), end(end_) {
  assert(start <= end);
}

inline bool Bin::operator==(const Bin &other) const noexcept {
  return this->chrom == other.chrom && this->start == other.start && this->end == other.end;
}
inline bool Bin::operator!=(const Bin &other) const noexcept { return !(*this == other); }

inline bool Bin::operator<(const Bin &other) const noexcept {
  if (this->chrom != other.chrom) {
    return this->chrom < other.chrom;
  }

  if (this->start != other.start) {
    return this->start < other.start;
  }

  return this->end < other.end;
}

inline bool Bin::operator<=(const Bin &other) const noexcept {
  if (this->chrom != other.chrom) {
    return this->chrom <= other.chrom;
  }

  if (this->start != other.start) {
    return this->start <= other.start;
  }

  return this->end <= other.end;
}

inline bool Bin::operator>(const Bin &other) const noexcept {
  if (this->chrom != other.chrom) {
    return this->chrom > other.chrom;
  }

  if (this->start != other.start) {
    return this->start > other.start;
  }

  return this->end > other.end;
}

inline bool Bin::operator>=(const Bin &other) const noexcept {
  if (this->chrom != other.chrom) {
    return this->chrom >= other.chrom;
  }

  if (this->start != other.start) {
    return this->start >= other.start;
  }

  return this->end >= other.end;
}

inline BinTableLazy::BinTableLazy(ChromosomeSet chroms, std::uint32_t bin_size)
    : _chroms(std::move(chroms)),
      _num_bins_prefix_sum(compute_num_bins_prefix_sum(_chroms, bin_size)),
      _bin_size(bin_size) {
  assert(bin_size != 0);
}

template <typename ChromIt>
inline BinTableLazy::BinTableLazy(ChromIt first_chrom, ChromIt last_chrom, std::uint32_t bin_size)
    : BinTableLazy(ChromosomeSet(first_chrom, last_chrom), bin_size) {}

template <typename ChromNameIt, typename ChromSizeIt>
inline BinTableLazy::BinTableLazy(ChromNameIt first_chrom_name, ChromNameIt last_chrom_name,
                                  ChromSizeIt first_chrom_size, std::uint32_t bin_size)
    : BinTableLazy(ChromosomeSet(first_chrom_name, last_chrom_name, first_chrom_size), bin_size) {}

inline std::size_t BinTableLazy::size() const noexcept {
  if (this->_num_bins_prefix_sum.empty()) {
    return 0;
  }
  return static_cast<std::size_t>(this->_num_bins_prefix_sum.back());
}

inline bool BinTableLazy::empty() const noexcept { return this->size() == 0; }

inline std::size_t BinTableLazy::num_chromosomes() const { return this->_chroms.size(); }

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

constexpr std::uint32_t BinTableLazy::iterator::bin_size() const noexcept {
  assert(this->_bin_table);
  return this->_bin_table->bin_size();
}

inline BinTable BinTableLazy::concretize() const {
  std::vector<const Chromosome *> chroms(this->size());
  std::vector<std::uint32_t> starts(this->size());
  std::vector<std::uint32_t> ends(this->size());

  std::size_t i = 0;
  for (const auto [chrom, start, end] : *this) {
    chroms[i] = &chrom;
    starts[i] = start;
    ends[i++] = end;
  }
  assert(i == chroms.size());

  return BinTable{chroms, starts, ends};
}

inline bool BinTableLazy::operator==(const BinTableLazy &other) const {
  return this->_bin_size == other._bin_size && this->_chroms == other._chroms;
}
inline bool BinTableLazy::operator!=(const BinTableLazy &other) const { return !(*this == other); }

inline BinTableLazy BinTableLazy::at(const Chromosome &chrom) const { return this->at(chrom.name); }

inline BinTableLazy BinTableLazy::at(std::string_view chrom_name) const {
  const auto match = this->_chroms.find(chrom_name);
  if (match == this->_chroms.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
  }

  return {ChromosomeSet{*match}, this->_bin_size};
}

inline BinTableLazy BinTableLazy::at(std::uint32_t chrom_id) const {
  auto match = this->_chroms.find(chrom_id);
  if (match == this->_chroms.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome with id {} not found"), chrom_id));
  }

  return {ChromosomeSet{*match}, this->_bin_size};
}

inline Bin BinTableLazy::bin_id_to_coords(std::uint64_t bin_id) const {
  auto match = std::upper_bound(this->_num_bins_prefix_sum.begin(),
                                this->_num_bins_prefix_sum.end(), bin_id);

  if (match == this->_num_bins_prefix_sum.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("bin id {} not found: out of range"), bin_id));
  }
  assert(match != this->_num_bins_prefix_sum.begin());
  --match;

  const auto chrom_id =
      static_cast<std::uint32_t>(std::distance(this->_num_bins_prefix_sum.begin(), match));
  assert(this->_chroms.contains(chrom_id));
  const auto &chrom = this->_chroms.at(chrom_id);

  const auto relative_bin_id = bin_id - *match;
  const auto start = static_cast<uint32_t>(relative_bin_id * this->bin_size());
  assert(start < chrom.size);
  const auto end = (std::min)(start + this->bin_size(), chrom.size);

  return {chrom, start, end};
}

inline std::uint64_t BinTableLazy::coord_to_bin_id(const Bin &bin) const {
  const auto match = this->_chroms.find(bin.chrom);
  if (match == this->_chroms.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), bin.chrom));
  }

  if (bin.end < bin.start) {
    throw std::logic_error(
        fmt::format(FMT_STRING("invalid coordinate: start > end: {} > {}"), bin.start, bin.end));
  }

  const auto chrom_id =
      conditional_static_cast<std::size_t>(std::distance(this->_chroms.begin(), match));
  const auto bin_offset = this->_num_bins_prefix_sum[chrom_id] - this->_num_bins_prefix_sum.front();

  return bin_offset + static_cast<std::uint64_t>(bin.start / this->bin_size());
}

inline std::uint64_t BinTableLazy::coord_to_bin_id(const Chromosome &chrom,
                                                   std::uint32_t pos) const {
  if (pos > chrom.size) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("position is greater than chromosome size: {} > {}"), pos, chrom));
  }
  return this->coord_to_bin_id(Bin{chrom, pos, (std::min)(pos + this->bin_size(), chrom.size)});
}

inline std::uint64_t BinTableLazy::coord_to_bin_id(std::string_view chrom_name,
                                                   std::uint32_t pos) const {
  auto match = this->_chroms.find(chrom_name);
  if (match == this->_chroms.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
  }

  if (pos > match->size) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("position is greater than chromosome size: {} > {}"), pos, *match));
  }

  return this->coord_to_bin_id(Bin{*match, pos, (std::min)(pos + this->bin_size(), match->size)});
}

inline std::uint64_t BinTableLazy::coord_to_bin_id(std::uint32_t chrom_id,
                                                   std::uint32_t pos) const {
  auto match = this->_chroms.find(chrom_id);
  if (match == this->_chroms.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome with id {} not found"), chrom_id));
  }

  if (pos > match->size) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("position is greater than chromosome size: {} > {}"), pos, *match));
  }

  return this->coord_to_bin_id(Bin{*match, pos, (std::min)(pos + this->bin_size(), match->size)});
}

inline std::vector<std::uint64_t> BinTableLazy::compute_num_bins_prefix_sum(
    const ChromosomeSet &chroms, std::uint32_t bin_size) {
  assert(bin_size != 0);

  std::vector<std::uint64_t> prefix_sum(chroms.size() + 1, 0);

  // I am using transform instead of inclusive_scan because the latter is not always available
  std::transform(chroms.begin(), chroms.end(), prefix_sum.begin() + 1,
                 [&, sum = std::uint64_t(0)](const Chromosome &chrom) mutable {
                   return sum += static_cast<std::uint64_t>((chrom.size + bin_size - 1) / bin_size);
                 });

  return prefix_sum;
}

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

inline auto BinTableLazy::iterator::operator*() const -> value_type {
  assert(this->_bin_table);

  const auto &chrom = this->chromosome();
  const auto bin_size = this->bin_size();

  const auto start = (std::min)(static_cast<std::uint32_t>(this->_idx) * bin_size, chrom.size);
  const auto end = (std::min)(start + bin_size, chrom.size);

  return value_type{chrom, start, end};
}

inline auto BinTableLazy::iterator::operator++() -> iterator & {
  assert(this->_bin_table);
  if (this->_chrom_id == std::numeric_limits<std::uint32_t>::max()) {
    return *this;
  }

  if (++this->_idx >= this->compute_num_bins()) {
    if (this->_chrom_id + 1 >= this->num_chromosomes()) {
      return *this = make_end_iterator(*this->_bin_table);
    }
    ++this->_chrom_id;
    this->_idx = 0;
  }

  return *this;
}

inline auto BinTableLazy::iterator::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

inline auto BinTableLazy::iterator::operator--() -> iterator & {
  assert(this->_bin_table);

  if (this->_idx == npos) {
    assert(*this == make_end_iterator(*this->_bin_table));
    this->_chrom_id = static_cast<std::uint32_t>(this->num_chromosomes() - 1);
    this->_idx = this->compute_num_bins() - 1;
    return *this;
  }

  if (this->_idx-- == 0) {
    this->_chrom_id--;
    this->_idx = this->compute_num_bins() - 1;
  }

  return *this;
}

inline auto BinTableLazy::iterator::operator--(int) -> iterator {
  auto it = *this;
  std::ignore = --(*this);
  return it;
}

inline const Chromosome &BinTableLazy::iterator::chromosome() const {
  return this->chromosome(this->_chrom_id);
}

inline const Chromosome &BinTableLazy::iterator::chromosome(std::uint32_t chrom_id) const {
  return this->_bin_table->chromosomes().at(chrom_id);
}

inline std::uint64_t BinTableLazy::iterator::compute_num_bins() const noexcept {
  assert(this->_bin_table);

  const auto chrom_size = this->chromosome().size;
  const auto bin_size = this->bin_size();

  return static_cast<std::uint64_t>((chrom_size + bin_size - 1) / bin_size);
}

inline std::size_t BinTableLazy::iterator::num_chromosomes() const noexcept {
  assert(this->_bin_table);

  return this->_bin_table->num_chromosomes();
}

}  // namespace coolerpp

constexpr auto fmt::formatter<coolerpp::Bin>::parse(format_parse_context &ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.end();
}

template <typename FormatContext>
inline auto fmt::formatter<coolerpp::Bin>::format(const coolerpp::Bin &b, FormatContext &ctx) const
    -> decltype(ctx.out()) {
  return fmt::format_to(ctx.out(), FMT_STRING("{}:{}-{}"), b.chrom.name, b.start, b.end);
}
