// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/index.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "coolerpp/bin_table.hpp"
#include "coolerpp/chromosome.hpp"

namespace coolerpp {

Index::Index(const BinTableLazy &bins, std::uint64_t nnz)
    : _bins(&bins), _idx(Index::init(bins.chromosomes(), _bins->bin_size())), _nnz(nnz) {
  assert(this->bin_size() != 0);
  _size = std::accumulate(_idx.begin(), _idx.end(), std::size_t(0),
                          [&](std::size_t sum, const auto &it) { return sum + it.size(); });
}

const ChromosomeSet &Index::chromosomes() const noexcept {
  assert(this->_bins);
  return this->_bins->chromosomes();
}

const BinTableLazy &Index::bins() const noexcept {
  assert(this->_bins);
  return *this->_bins;
}

std::size_t Index::num_chromosomes() const noexcept {
  assert(this->_idx.size() == this->_bins->num_chromosomes());
  return this->_idx.size();
}

std::size_t Index::size(std::string_view chrom_name) const {
  const auto chrom_id = this->chromosomes().get_id(chrom_name);
  return this->size(chrom_id);
}

std::size_t Index::size(std::uint32_t chrom_id) const {
  this->validate_chrom_id(chrom_id);
  return this->at(chrom_id).size();
}

std::uint32_t Index::bin_size() const noexcept {
  assert(this->_bins);
  return this->_bins->bin_size();
}

auto Index::begin() const noexcept -> const_iterator { return iterator{this}; }
auto Index::end() const noexcept -> const_iterator { return iterator::make_end_iterator(this); }

auto Index::cbegin() const noexcept -> const_iterator { return this->begin(); }
auto Index::cend() const noexcept -> const_iterator { return this->end(); }

auto Index::at(std::string_view chrom_name) const -> const mapped_type & {
  const auto chrom_id = this->chromosomes().get_id(chrom_name);
  return this->_idx.at(chrom_id);
}

auto Index::at(std::uint32_t chrom_id) -> mapped_type & {
  this->validate_chrom_id(chrom_id);
  return this->_idx.at(chrom_id);
}

auto Index::at(std::string_view chrom_name) -> mapped_type & {
  const auto chrom_id = this->chromosomes().get_id(chrom_name);
  return this->_idx.at(chrom_id);
}

auto Index::at(std::uint32_t chrom_id) const -> const mapped_type & {
  this->validate_chrom_id(chrom_id);
  return this->_idx.at(chrom_id);
}

std::uint64_t Index::get_offset_by_bin_id(std::uint64_t bin_id) const {
  if (bin_id == this->size()) {
    return this->_idx.back().back();
  }
  const auto &coords = this->_bins->bin_id_to_coords(bin_id);
  return this->get_offset_by_pos(coords.chrom, coords.bin_start);
}

std::uint64_t Index::get_offset_by_pos(const Chromosome &chrom, std::uint32_t pos) const {
  return this->get_offset_by_pos(chrom.name, pos);
}

std::uint64_t Index::get_offset_by_pos(std::string_view chrom_name, std::uint32_t pos) const {
  const auto row_idx = pos / this->bin_size();
  return this->get_offset_by_row_idx(this->chromosomes().get_id(chrom_name), row_idx);
}

std::uint64_t Index::get_offset_by_pos(std::uint32_t chrom_id, std::uint32_t pos) const {
  const auto row_idx = pos / this->bin_size();
  return this->get_offset_by_row_idx(chrom_id, row_idx);
}

std::uint64_t Index::get_offset_by_row_idx(std::uint32_t chrom_id, std::size_t row_idx) const {
  const auto &offsets = this->at(chrom_id);
  if (row_idx >= offsets.size()) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("invalid row_index {}: row maps outside of chromosome {}"), row_idx,
                    this->chromosomes().at(chrom_id)));
  }
  return offsets[row_idx];
}

void Index::set_offset_by_bin_id(std::uint64_t bin_id, std::uint64_t offset) {
  const auto &coords = this->_bins->bin_id_to_coords(bin_id);
  this->set_offset_by_pos(coords.chrom, coords.bin_start, offset);
}

void Index::set_offset_by_pos(const Chromosome &chrom, std::uint32_t pos, std::uint64_t offset) {
  this->set_offset_by_pos(chrom.name, pos, offset);
}

void Index::set_offset_by_pos(std::string_view chrom_name, std::uint32_t pos,
                              std::uint64_t offset) {
  const auto row_idx = pos / this->bin_size();
  this->set_offset_by_row_idx(this->chromosomes().get_id(chrom_name), row_idx, offset);
}

void Index::set_offset_by_pos(std::uint32_t chrom_id, std::uint32_t pos, std::uint64_t offset) {
  const auto row_idx = pos / this->bin_size();
  this->set_offset_by_row_idx(chrom_id, row_idx, offset);
}

void Index::set_offset_by_row_idx(std::uint32_t chrom_id, std::size_t row_idx,
                                  std::uint64_t offset) {
  auto &offsets = this->at(chrom_id);
  assert(row_idx < offsets.size());
  offsets[row_idx] = offset;
}

void Index::validate() const {
  std::for_each(this->chromosomes().begin(), this->chromosomes().end(),
                [this](const Chromosome &chrom) { this->validate(chrom); });
}

std::uint64_t &Index::nnz() noexcept { return this->_nnz; }

std::vector<std::uint64_t> Index::compute_chrom_offsets() const {
  std::vector<std::uint64_t> buff(this->num_chromosomes());
  this->compute_chrom_offsets(buff);
  return buff;
}

std::uint64_t Index::get_bin1_offset(std::string_view chrom_name) const {
  return this->at(chrom_name).front();
}

std::uint64_t Index::get_bin1_offset(std::uint32_t chrom_id) const {
  return this->at(chrom_id).front();
}

void Index::finalize(std::uint64_t nnz) {
  this->_nnz = nnz;
  auto fill_value = nnz;

  std::for_each(this->_idx.rbegin(), this->_idx.rend(), [&](OffsetVect &offsets) {
    std::transform(offsets.rbegin(), offsets.rend(), offsets.rbegin(), [&fill_value](auto &offset) {
      if (offset == Index::offset_not_set_value) {
        return fill_value;
      }

      return fill_value = offset;
    });
  });
  assert(this->_idx[0][0] == 0 || this->_idx[0][0] == this->_idx[0][1]);
  this->_idx[0][0] = 0;
}

void Index::compute_chrom_offsets(std::vector<std::uint64_t> &buff) const noexcept {
  buff.resize(this->num_chromosomes() + 1);
  buff[0] = 0;

  std::transform(this->_idx.begin(), this->_idx.end(), buff.begin() + 1,
                 [offset = std::uint64_t(0)](const auto &it) mutable {
                   return offset += conditional_static_cast<std::uint64_t>(it.size());
                 });
}

void Index::validate_chrom_id(std::uint32_t chrom_id) const {
  if (static_cast<std::size_t>(chrom_id) >= this->num_chromosomes()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome with id {} not found"), chrom_id));
  }
}

auto Index::init(const ChromosomeSet &chroms, std::uint32_t bin_size) -> MapT {
  assert(!chroms.empty());
  assert(bin_size != 0);

  MapT idx(chroms.size());
  std::transform(chroms.begin(), chroms.end(), idx.begin(), [&](const Chromosome &chrom) {
    const auto num_bins = (chrom.size + bin_size - 1) / bin_size;
    return std::vector<std::uint64_t>(num_bins, Index::offset_not_set_value);
  });

  return idx;
}

void Index::validate(const Chromosome &chrom) const {
  try {
    const auto chrom_id = this->chromosomes().get_id(chrom);
    const auto &offsets = this->at(chrom_id);
    if (chrom_id == 0) {
      if (offsets.front() != 0) {
        throw std::runtime_error("first offset is not zero");
      }
    } else {
      const auto &prev_offsets = this->at(chrom_id - 1);
      if (offsets.front() < prev_offsets.back()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("offsets are not in ascending order: offset for "
                                   "bin {}:{}-{} should be >= {}, found {}"),
                        chrom.name, 0, this->bin_size(), prev_offsets.back(), offsets.front()));
      }
    }

    if (const auto it = std::is_sorted_until(offsets.begin(), offsets.end()); it != offsets.end()) {
      const auto i = std::distance(offsets.begin(), it);
      assert(i != 0);
      throw std::runtime_error(
          fmt::format(FMT_STRING("offsets are not in ascending order: pixels/bin1_offset[{}]={} > "
                                 "pixels/bin1_offset[{}]={}\n"),
                      i - 1, *(it - 1), i, *(it)));
    }

    if (this->_nnz != 0) {
      auto match = std::find_if(offsets.begin(), offsets.end(),
                                [this](const auto offset) { return offset > this->_nnz; });
      if (match != offsets.end()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("invalid offset {}: offset is greater than nnz ({} > {})"),
                        *match, *match, this->_nnz));
      }
    }

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("{} index is corrupted or incomplete: {}"), chrom.name, e.what()));
  }
}

// NOLINTNEXTLINE
Index::iterator::iterator(const Index *idx) : _idx(idx), _chrom_id(0), _offset_idx(0) {
  assert(idx);
}

bool Index::iterator::operator==(const iterator &other) const noexcept {
  return this->_idx == other._idx && this->_chrom_id == other._chrom_id &&
         this->_offset_idx == other._offset_idx;
}

bool Index::iterator::operator!=(const iterator &other) const noexcept { return !(*this == other); }

auto Index::iterator::operator*() const -> value_type {
  assert(this->_idx);
  if (this->_chrom_id > this->last_chrom_id()) {
    return this->_idx->_nnz;
  }
  return this->get_offsets()[this->_offset_idx];
}

auto Index::iterator::operator++() -> iterator & {
  if (this->_chrom_id > this->last_chrom_id()) {
    return *this = make_end_iterator(this->_idx);
  }

  if (++this->_offset_idx >= this->get_offsets().size()) {
    if (++this->_chrom_id > this->last_chrom_id()) {
      return *this;  // Next dereference returns the index size
    }

    this->_offset_idx = 0;
  }

  return *this;
}

auto Index::iterator::operator++(int) -> iterator {
  auto it = *this;
  ++(*this);
  return it;
}

auto Index::iterator::make_end_iterator(const Index *idx) -> iterator {
  assert(idx);

  iterator it{};

  it._idx = idx;
  it._chrom_id = it.last_chrom_id() + 1;
  it._offset_idx = npos;

  return it;
}

std::uint32_t Index::iterator::last_chrom_id() const noexcept {
  assert(this->_idx);
  if (this->_idx->num_chromosomes() == 0) {
    return 0;
  }

  return static_cast<std::uint32_t>(this->_idx->num_chromosomes() - 1);
}

auto Index::iterator::get_offsets() const noexcept -> const OffsetVect & {
  assert(this->_chrom_id < static_cast<std::uint32_t>(this->_idx->size()));
  return this->_idx->_idx[static_cast<std::size_t>(this->_chrom_id)];
}

}  // namespace coolerpp
