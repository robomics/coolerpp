// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <utility>

#include "coolerpp/bin_table.hpp"
#include "coolerpp/index.hpp"
#include "coolerpp/internal/numeric_utils.hpp"

namespace coolerpp {

template <class N>
inline PixelSelector<N>::PixelSelector(const Index &index, const Dataset &pixels_bin1_id,
                                       const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                                       std::shared_ptr<PixelCoordinates> coords) noexcept
    : PixelSelector(index, pixels_bin1_id, pixels_bin2_id, pixels_count, coords,
                    std::move(coords)) {}

template <class N>
inline PixelSelector<N>::PixelSelector(const Index &index, const Dataset &pixels_bin1_id,
                                       const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                                       std::shared_ptr<PixelCoordinates> coord1,
                                       std::shared_ptr<PixelCoordinates> coord2) noexcept
    : _coord1(std::move(coord1)),
      _coord2(std::move(coord2)),
      _index(&index),
      _pixels_bin1_id(&pixels_bin1_id),
      _pixels_bin2_id(&pixels_bin2_id),
      _pixels_count(&pixels_count) {}

template <class N>
inline PixelSelector<N>::PixelSelector(const Index &index, const Dataset &pixels_bin1_id,
                                       const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                                       PixelCoordinates coords) noexcept
    : PixelSelector<N>(index, pixels_bin1_id, pixels_bin2_id, pixels_count,
                       std::make_shared<PixelCoordinates>(std::move(coords))) {}

template <class N>
inline PixelSelector<N>::PixelSelector(const Index &index, const Dataset &pixels_bin1_id,
                                       const Dataset &pixels_bin2_id,
                                       const Dataset &pixels_count) noexcept
    : PixelSelector<N>(index, pixels_bin1_id, pixels_bin2_id, pixels_count, nullptr, nullptr) {}

template <class N>
inline PixelSelector<N>::PixelSelector(const Index &index, const Dataset &pixels_bin1_id,
                                       const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                                       PixelCoordinates coord1, PixelCoordinates coord2) noexcept
    : PixelSelector(index, pixels_bin1_id, pixels_bin2_id, pixels_count,
                    std::make_shared<PixelCoordinates>(std::move(coord1)),
                    std::make_shared<PixelCoordinates>(std::move(coord2))) {}

template <class N>
inline auto PixelSelector<N>::begin() const -> iterator {
  return this->cbegin();
}

template <class N>
inline auto PixelSelector<N>::end() const -> iterator {
  return this->cend();
}

template <class N>
inline auto PixelSelector<N>::cbegin() const -> iterator {
  if (!this->_coord1) {
    assert(!this->_coord2);
    return iterator{*this->_index, *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                    *this->_pixels_count};
  }

  return iterator{*this->_index,        *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                  *this->_pixels_count, this->_coord1,          this->_coord2};
}

template <class N>
inline auto PixelSelector<N>::cend() const -> iterator {
  return iterator::at_end(*this->_index, *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                          *this->_pixels_count, this->_coord1, this->_coord2);
}

template <class N>
constexpr const PixelCoordinates &PixelSelector<N>::coord1() const noexcept {
  return this->_coord1;
}

template <class N>
constexpr const PixelCoordinates &PixelSelector<N>::coord2() const noexcept {
  return this->_coord2;
}

template <class N>
inline PixelCoordinates PixelSelector<N>::parse_query(const BinTableLazy &bins,
                                                      std::string_view query) {
  if (bins.chromosomes().contains(query)) {
    const auto &chrom = bins.chromosomes().at(query);
    return {bins, chrom, 0, chrom.size};
  }

  const auto p1 = query.find_last_of(':');
  const auto p2 = query.find_last_of('-');

  if (p1 == std::string_view::npos || p2 == std::string_view::npos || p1 > p2) {
    throw std::runtime_error(fmt::format(FMT_STRING("query \"{}\" is malformed"), query));
  }

  const auto chrom_name = query.substr(0, p1);
  const auto start_pos_str = query.substr(p1 + 1, p2 - (p1 + 1));
  const auto end_pos_str = query.substr(p2 + 1);

  if (!bins.chromosomes().contains(chrom_name)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid chromosome \"{}\" in query \"{}\""), chrom_name, query));
  }

  if (start_pos_str.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("query \"{}\" is malformed: missing start position"), query));
  }

  if (end_pos_str.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("query \"{}\" is malformed: missing end position"), query));
  }

  const auto &chrom = bins.chromosomes().at(chrom_name);
  std::uint32_t start_pos{};
  std::uint32_t end_pos{};

  try {
    internal::parse_numeric_or_throw(start_pos_str, start_pos);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid start position \"{}\" in query \"{}\": {}"), start_pos_str,
                    query, e.what()));
  }
  try {
    internal::parse_numeric_or_throw(end_pos_str, end_pos);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid end position \"{}\" in query \"{}\": {}"), end_pos_str,
                    query, e.what()));
  }

  if (end_pos > chrom.size) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid end position \"{}\" in query \"{}\": end position is "
                               "greater than the chromosome size ({} > {})"),
                    end_pos, query, end_pos, chrom.size));
  }

  if (start_pos >= end_pos) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid query \"{}\": query end position should be "
                               "greater than the start position ({} >= {})"),
                    query, start_pos, end_pos));
  }

  end_pos -= std::min(end_pos, 1U);
  return {bins, chrom, start_pos, end_pos};
}

template <class N>
inline PixelSelector<N>::iterator::iterator(const Index &index, const Dataset &pixels_bin1_id,
                                            const Dataset &pixels_bin2_id,
                                            const Dataset &pixels_count)
    : _index(&index),
      _coord1(nullptr),
      _coord2(nullptr),
      _bin1_id_it(pixels_bin1_id.begin<std::uint64_t>()),
      _bin2_id_it(pixels_bin2_id.begin<std::uint64_t>()),
      _bin2_id_last(pixels_bin2_id.end<std::uint64_t>()),
      _count_it(pixels_count.begin<N>()) {}

template <class N>
inline PixelSelector<N>::iterator::iterator(const Index &index, const Dataset &pixels_bin1_id,
                                            const Dataset &pixels_bin2_id,
                                            const Dataset &pixels_count,
                                            std::shared_ptr<PixelCoordinates> coord1,
                                            std::shared_ptr<PixelCoordinates> coord2)
    : _index(&index), _coord1(std::move(coord1)), _coord2(std::move(coord2)) {
  assert(_coord1);
  assert(_coord2);
  assert(_coord1->bin1_id() <= _coord1->bin2_id());
  assert(_coord2->bin1_id() <= _coord2->bin2_id());

  auto offset = _index->get_offset_by_bin_id(_coord1->bin1_id());
  _bin1_id_it = pixels_bin1_id.make_iterator_at_offset<std::uint64_t>(offset);
  _bin2_id_it = pixels_bin2_id.make_iterator_at_offset<std::uint64_t>(offset);
  _bin2_id_last = pixels_bin2_id.end<std::uint64_t>();
  _count_it = pixels_count.make_iterator_at_offset<N>(offset);

  this->jump(_coord1->bin1_id(), _coord2->bin1_id());

  auto it = iterator::at_end(index, pixels_bin1_id, pixels_bin2_id, pixels_count, _coord1, _coord2);
  _bin2_id_last = it._bin2_id_it;
  assert(_bin2_id_it <= _bin2_id_last);
}

template <class N>
inline auto PixelSelector<N>::iterator::at_end(const Index &index, const Dataset &pixels_bin1_id,
                                               const Dataset &pixels_bin2_id,
                                               const Dataset &pixels_count) -> iterator {
  iterator it{};
  it._index = &index;
  it._bin1_id_it = pixels_bin1_id.end<std::uint64_t>();
  it._bin2_id_it = pixels_bin2_id.end<std::uint64_t>();
  it._bin2_id_last = it._bin2_id_it;
  it._count_it = pixels_count.end<N>();

  return it;
}

template <class N>
inline auto PixelSelector<N>::iterator::at_end(const Index &index, const Dataset &pixels_bin1_id,
                                               const Dataset &pixels_bin2_id,
                                               const Dataset &pixels_count,
                                               std::shared_ptr<PixelCoordinates> coord1,
                                               std::shared_ptr<PixelCoordinates> coord2)
    -> iterator {
  if (!coord1 && !coord2) {
    return at_end(index, pixels_bin1_id, pixels_bin2_id, pixels_count);
  }
  assert(!!coord1);
  assert(!!coord2);

  iterator it{};

  it._index = &index;
  it._coord1 = std::move(coord1);
  it._coord2 = std::move(coord2);

  const auto offset = index.get_offset_by_bin_id(it._coord1->bin2_id());
  it._bin1_id_it = pixels_bin1_id.make_iterator_at_offset<std::uint64_t>(offset);
  it._bin2_id_it = pixels_bin2_id.make_iterator_at_offset<std::uint64_t>(offset);
  it._bin2_id_last = pixels_bin2_id.end<std::uint64_t>();
  it._count_it = pixels_count.make_iterator_at_offset<N>(offset);

  it.jump(it._coord1->bin2_id(), std::max(it._coord1->bin2_id(), it._coord2->bin2_id()));

  it._bin2_id_last = ++it._bin2_id_it;

  return it;
}

template <class N>
constexpr bool PixelSelector<N>::iterator::operator==(const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it == other._bin2_id_it;
}
template <class N>
constexpr bool PixelSelector<N>::iterator::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <class N>
constexpr bool PixelSelector<N>::iterator::operator<(const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it < other._bin2_id_it;
}
template <class N>
constexpr bool PixelSelector<N>::iterator::operator<=(const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it <= other._bin2_id_it;
}

template <class N>
constexpr bool PixelSelector<N>::iterator::operator>(const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it > other._bin2_id_it;
}
template <class N>
constexpr bool PixelSelector<N>::iterator::operator>=(const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it >= other._bin2_id_it;
}

template <class N>
inline auto PixelSelector<N>::iterator::operator*() const -> value_type {
  assert(this->_index);
  assert(this->_bin2_id_it < this->_bin2_id_last);
  // clang-format off
  return {PixelCoordinates{this->_index->bins(),
                           *this->_bin1_id_it,
                           *this->_bin2_id_it},
          *this->_count_it};
  // clang-format on
}

template <class N>
inline auto PixelSelector<N>::iterator::operator++() -> iterator & {
  if (!this->_coord1) {
    assert(!this->_coord2);
    std::ignore = ++this->_bin1_id_it;
    std::ignore = ++this->_bin2_id_it;
    std::ignore = ++this->_count_it;
    return *this;
  }

  if (this->_bin2_id_it < this->_bin2_id_last && *this->_bin2_id_it >= this->_coord2->bin2_id()) {
    const auto row = *this->_bin1_id_it + 1;
    if (row <= this->_coord1->bin2_id()) {
      const auto col = std::max(row, this->_coord2->bin1_id());
      this->jump(row, col);
      return *this;
    }
  }

  std::ignore = ++this->_bin1_id_it;
  std::ignore = ++this->_bin2_id_it;
  std::ignore = ++this->_count_it;

  return *this;
}

template <class N>
inline auto PixelSelector<N>::iterator::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <class N>
inline void PixelSelector<N>::iterator::jump_to_row(std::uint64_t bin_id) {
  assert(this->_index);
  assert(bin_id <= this->_index->bins().size());

  auto offset = this->_index->get_offset_by_bin_id(bin_id);

  const auto &bin2_id_dset = _bin2_id_it.dataset();
  const auto &bin1_id_dset = _bin1_id_it.dataset();
  const auto &count_dset = _count_it.dataset();

  this->_bin1_id_it = bin1_id_dset.make_iterator_at_offset<std::uint64_t>(offset);
  this->_bin2_id_it = bin2_id_dset.make_iterator_at_offset<std::uint64_t>(offset);
  this->_count_it = count_dset.template make_iterator_at_offset<N>(offset);

  assert(this->_bin2_id_it <= this->_bin2_id_last);
}

template <class N>
inline void PixelSelector<N>::iterator::jump_to_col(std::uint64_t bin_id) {
  assert(this->_index);
  assert(bin_id <= this->_index->bins().size());

  const auto current_row = *this->_bin1_id_it;
  const auto next_row = current_row + 1;

  const auto offset1 = this->_index->get_offset_by_bin_id(current_row);
  const auto offset2 = [&]() {
    const auto offset = this->_index->get_offset_by_bin_id(next_row);
    if (offset == 0) {
      return offset;
    }
    return std::max(offset1, offset - 1);
  }();

  assert(offset1 <= offset2);
  if (offset1 == offset2) {
    return;
  }

  const auto &bin2_id_dset = this->_bin2_id_it.dataset();
  this->_bin2_id_it = std::lower_bound(
      bin2_id_dset.template make_iterator_at_offset<std::uint64_t>(offset1),
      bin2_id_dset.template make_iterator_at_offset<std::uint64_t>(offset2), bin_id);

  const auto &bin1_id_dset = this->_bin1_id_it.dataset();
  const auto &count_dset = this->_count_it.dataset();
  const auto offset = this->_bin2_id_it.h5_offset();

  this->_bin1_id_it = bin1_id_dset.template make_iterator_at_offset<std::uint64_t>(offset);
  this->_count_it = count_dset.template make_iterator_at_offset<N>(offset);

  assert(*this->_bin1_id_it == current_row);
  assert(this->_bin2_id_it <= this->_bin2_id_last);
}

template <class N>
inline void PixelSelector<N>::iterator::jump(std::uint64_t bin1_id, std::uint64_t bin2_id) {
  assert(bin1_id <= bin2_id);

  this->jump_to_row(bin1_id);
  if (bin2_id != bin1_id) {
    this->jump_to_col(bin2_id);
  }
}

}  // namespace coolerpp
