// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/bin_table.hpp"

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

namespace coolerpp {

bool Bin::operator==(const Bin &other) const noexcept {
  return this->chrom == other.chrom && this->bin_start == other.bin_start &&
         this->bin_end == other.bin_end;
}
bool Bin::operator!=(const Bin &other) const noexcept { return !(*this == other); }

bool Bin::operator<(const Bin &other) const noexcept {
  if (this->chrom != other.chrom) {
    return this->chrom < other.chrom;
  }

  if (this->bin_start != other.bin_start) {
    return this->bin_start < other.bin_start;
  }

  return this->bin_end < other.bin_end;
}

bool Bin::operator<=(const Bin &other) const noexcept {
  if (this->chrom != other.chrom) {
    return this->chrom <= other.chrom;
  }

  if (this->bin_start != other.bin_start) {
    return this->bin_start <= other.bin_start;
  }

  return this->bin_end <= other.bin_end;
}

bool Bin::operator>(const Bin &other) const noexcept {
  if (this->chrom != other.chrom) {
    return this->chrom > other.chrom;
  }

  if (this->bin_start != other.bin_start) {
    return this->bin_start > other.bin_start;
  }

  return this->bin_end > other.bin_end;
}

bool Bin::operator>=(const Bin &other) const noexcept {
  if (this->chrom != other.chrom) {
    return this->chrom >= other.chrom;
  }

  if (this->bin_start != other.bin_start) {
    return this->bin_start >= other.bin_start;
  }

  return this->bin_end >= other.bin_end;
}

BinTableLazy::BinTableLazy(ChromosomeSet chroms, std::uint32_t bin_size)
    : _chroms(std::move(chroms)),
      _num_bins_prefix_sum(compute_num_bins_prefix_sum(_chroms, bin_size)),
      _bin_size(bin_size) {
  assert(bin_size != 0);
}

std::size_t BinTableLazy::size() const noexcept {
  if (this->_num_bins_prefix_sum.empty()) {
    return 0;
  }
  return static_cast<std::size_t>(this->_num_bins_prefix_sum.back());
}

bool BinTableLazy::empty() const noexcept { return this->size() == 0; }

std::size_t BinTableLazy::num_chromosomes() const { return this->_chroms.size(); }

BinTable BinTableLazy::concretize() const {
  std::vector<const Chromosome *> chroms(this->size());
  std::vector<std::uint32_t> bin_starts(this->size());
  std::vector<std::uint32_t> bin_ends(this->size());

  std::size_t i = 0;
  for (const auto [chrom, bin_start, bin_end] : *this) {
    chroms[i] = &chrom;
    bin_starts[i] = bin_start;
    bin_ends[i++] = bin_end;
  }
  assert(i == chroms.size());

  return BinTable{chroms, bin_starts, bin_ends};
}

bool BinTableLazy::operator==(const BinTableLazy &other) const {
  return this->_bin_size == other._bin_size && this->_chroms == other._chroms;
}
bool BinTableLazy::operator!=(const BinTableLazy &other) const { return !(*this == other); }

BinTableLazy BinTableLazy::at(const Chromosome &chrom) const { return this->at(chrom.name); }

BinTableLazy BinTableLazy::at(std::string_view chrom_name) const {
  const auto match = this->_chroms.find(chrom_name);
  if (match == this->_chroms.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
  }

  return {ChromosomeSet{*match}, this->_bin_size};
}

BinTableLazy BinTableLazy::at(std::uint32_t chrom_id) const {
  auto match = this->_chroms.find(chrom_id);
  if (match == this->_chroms.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome with id {} not found"), chrom_id));
  }

  return {ChromosomeSet{*match}, this->_bin_size};
}

Bin BinTableLazy::bin_id_to_coords(std::uint64_t bin_id) const {
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
  const auto bin_start = static_cast<uint32_t>(relative_bin_id * this->bin_size());
  assert(bin_start < chrom.size);
  const auto bin_end = std::min(bin_start + this->bin_size(), chrom.size);

  return {chrom, bin_start, bin_end};
}

std::uint64_t BinTableLazy::coord_to_bin_id(const Bin &bin) const {
  const auto match = this->_chroms.find(bin.chrom);
  if (match == this->_chroms.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), bin.chrom));
  }

  if (bin.bin_end < bin.bin_start) {
    throw std::logic_error(
        fmt::format(FMT_STRING("invalid coordinate: bin_start > bin_end: {} > {}"), bin.bin_start,
                    bin.bin_end));
  }

  const auto chrom_id =
      conditional_static_cast<std::size_t>(std::distance(this->_chroms.begin(), match));
  const auto bin_offset = this->_num_bins_prefix_sum[chrom_id] - this->_num_bins_prefix_sum.front();

  return bin_offset + static_cast<std::uint64_t>(bin.bin_start / this->bin_size());
}

std::uint64_t BinTableLazy::coord_to_bin_id(const Chromosome &chrom, std::uint32_t pos) const {
  if (pos > chrom.size) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("position is greater than chromosome size: {} > {}"), pos, chrom));
  }
  return this->coord_to_bin_id(Bin{chrom, pos, std::min(pos + this->bin_size(), chrom.size)});
}

std::uint64_t BinTableLazy::coord_to_bin_id(std::string_view chrom_name, std::uint32_t pos) const {
  auto match = this->_chroms.find(chrom_name);
  if (match == this->_chroms.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
  }

  if (pos > match->size) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("position is greater than chromosome size: {} > {}"), pos, *match));
  }

  return this->coord_to_bin_id(Bin{*match, pos, std::min(pos + this->bin_size(), match->size)});
}

std::uint64_t BinTableLazy::coord_to_bin_id(std::uint32_t chrom_id, std::uint32_t pos) const {
  auto match = this->_chroms.find(chrom_id);
  if (match == this->_chroms.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome with id {} not found"), chrom_id));
  }

  if (pos > match->size) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("position is greater than chromosome size: {} > {}"), pos, *match));
  }

  return this->coord_to_bin_id(Bin{*match, pos, std::min(pos + this->bin_size(), match->size)});
}

std::vector<std::uint64_t> BinTableLazy::compute_num_bins_prefix_sum(const ChromosomeSet &chroms,
                                                                     std::uint32_t bin_size) {
  assert(bin_size != 0);

  std::vector<std::uint64_t> prefix_sum(chroms.size() + 1, 0);

  // I am using transform instead of inclusive_scan because the latter is not always available
  std::transform(chroms.begin(), chroms.end(), prefix_sum.begin() + 1,
                 [&, sum = std::uint64_t(0)](const Chromosome &chrom) mutable {
                   return sum += static_cast<std::uint64_t>((chrom.size + bin_size - 1) / bin_size);
                 });

  return prefix_sum;
}

auto BinTableLazy::iterator::operator*() const -> value_type {
  assert(this->_bin_table);

  const auto &chrom = this->chromosome();
  const auto bin_size = this->bin_size();

  const auto start = std::min(static_cast<std::uint32_t>(this->_idx) * bin_size, chrom.size);
  const auto end = std::min(start + bin_size, chrom.size);

  return value_type{chrom, start, end};
}

auto BinTableLazy::iterator::operator++() -> iterator & {
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

auto BinTableLazy::iterator::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

auto BinTableLazy::iterator::operator--() -> iterator & {
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

auto BinTableLazy::iterator::operator--(int) -> iterator {
  auto it = *this;
  std::ignore = --(*this);
  return it;
}

const Chromosome &BinTableLazy::iterator::chromosome() const {
  return this->chromosome(this->_chrom_id);
}

const Chromosome &BinTableLazy::iterator::chromosome(std::uint32_t chrom_id) const {
  return this->_bin_table->chromosomes().at(chrom_id);
}

std::uint64_t BinTableLazy::iterator::compute_num_bins() const noexcept {
  assert(this->_bin_table);

  const auto chrom_size = this->chromosome().size;
  const auto bin_size = this->bin_size();

  return static_cast<std::uint64_t>((chrom_size + bin_size - 1) / bin_size);
}

std::size_t BinTableLazy::iterator::num_chromosomes() const noexcept {
  assert(this->_bin_table);

  return this->_bin_table->num_chromosomes();
}

}  // namespace coolerpp
