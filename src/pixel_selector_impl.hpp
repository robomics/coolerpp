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
inline PixelSelector<N>::PixelSelector(std::shared_ptr<const Index> index,
                                       const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                                       const Dataset &pixels_count,
                                       const std::shared_ptr<PixelCoordinates> &coords) noexcept
    : PixelSelector(std::move(index), pixels_bin1_id, pixels_bin2_id, pixels_count, coords,
                    coords) {}

template <class N>
inline PixelSelector<N>::PixelSelector(std::shared_ptr<const Index> index,
                                       const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                                       const Dataset &pixels_count,
                                       std::shared_ptr<PixelCoordinates> coord1,
                                       std::shared_ptr<PixelCoordinates> coord2) noexcept
    : _coord1(std::move(coord1)),
      _coord2(std::move(coord2)),
      _index(std::move(index)),
      _pixels_bin1_id(&pixels_bin1_id),
      _pixels_bin2_id(&pixels_bin2_id),
      _pixels_count(&pixels_count) {
  assert(_index);
}

template <class N>
inline PixelSelector<N>::PixelSelector(std::shared_ptr<const Index> index,
                                       const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                                       const Dataset &pixels_count,
                                       PixelCoordinates coords) noexcept
    : PixelSelector<N>(std::move(index), pixels_bin1_id, pixels_bin2_id, pixels_count,
                       std::make_shared<PixelCoordinates>(std::move(coords))) {}

template <class N>
inline PixelSelector<N>::PixelSelector(std::shared_ptr<const Index> index,
                                       const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                                       const Dataset &pixels_count) noexcept
    : PixelSelector<N>(std::move(index), pixels_bin1_id, pixels_bin2_id, pixels_count, nullptr,
                       nullptr) {}

template <class N>
inline PixelSelector<N>::PixelSelector(std::shared_ptr<const Index> index,
                                       const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                                       const Dataset &pixels_count, PixelCoordinates coord1,
                                       PixelCoordinates coord2) noexcept
    : PixelSelector(std::move(index), pixels_bin1_id, pixels_bin2_id, pixels_count,
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
    return iterator{this->_index, *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                    *this->_pixels_count};
  }

  return iterator{this->_index,         *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                  *this->_pixels_count, this->_coord1,          this->_coord2};
}

template <class N>
inline auto PixelSelector<N>::cend() const -> iterator {
  return iterator::at_end(this->_index, *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                          *this->_pixels_count, this->_coord1, this->_coord2);
}

template <class N>
constexpr const PixelCoordinates &PixelSelector<N>::coord1() const noexcept {
  assert(this->_coord1);
  return *this->_coord1;
}

template <class N>
constexpr const PixelCoordinates &PixelSelector<N>::coord2() const noexcept {
  assert(this->_coord1);
  return *this->_coord2;
}

template <class N>
inline PixelCoordinates PixelSelector<N>::parse_query(std::shared_ptr<const BinTableLazy> bins,
                                                      std::string_view query) {
  assert(bins);
  const auto &chroms = bins->chromosomes();
  if (chroms.contains(query)) {
    const auto &chrom = chroms.at(query);
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

  if (!chroms.contains(chrom_name)) {
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

  const auto &chrom = chroms.at(chrom_name);
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
inline PixelSelector<N>::iterator::iterator(std::shared_ptr<const Index> index,
                                            const Dataset &pixels_bin1_id,
                                            const Dataset &pixels_bin2_id,
                                            const Dataset &pixels_count)
    : _index(std::move(index)),
      _coord1(nullptr),
      _coord2(nullptr),
      _bin1_id_it(pixels_bin1_id.begin<std::uint64_t>()),
      _bin2_id_it(pixels_bin2_id.begin<std::uint64_t>()),
      _bin2_id_last(pixels_bin2_id.end<std::uint64_t>()),
      _count_it(pixels_count.begin<N>()) {}

template <class N>
inline PixelSelector<N>::iterator::iterator(std::shared_ptr<const Index> index,
                                            const Dataset &pixels_bin1_id,
                                            const Dataset &pixels_bin2_id,
                                            const Dataset &pixels_count,
                                            std::shared_ptr<PixelCoordinates> coord1,
                                            std::shared_ptr<PixelCoordinates> coord2)
    : _index(std::move(index)), _coord1(std::move(coord1)), _coord2(std::move(coord2)) {
  assert(_coord1);
  assert(_coord2);
  assert(_coord1->bin1_id() <= _coord1->bin2_id());
  assert(_coord2->bin1_id() <= _coord2->bin2_id());

  // Set iterator to the first row overlapping the query (i.e. the first bin overlapping coord1)
  auto offset = _index->get_offset_by_bin_id(_coord1->bin1_id());
  _bin1_id_it = pixels_bin1_id.make_iterator_at_offset<std::uint64_t>(offset);
  _bin2_id_it = pixels_bin2_id.make_iterator_at_offset<std::uint64_t>(offset);
  _count_it = pixels_count.make_iterator_at_offset<N>(offset);

  // Init last iterator. at_end will return an iterator pointing to the pixel past the last pixel
  // overlapping the query
  auto it = iterator::at_end(this->_index, pixels_bin1_id, pixels_bin2_id, pixels_count, _coord1,
                             _coord2);
  _bin2_id_last = it._bin2_id_last;
  assert(_bin2_id_it <= _bin2_id_last);

  if (_bin2_id_it.h5_offset() == pixels_bin2_id.size()) {
    return;
  }

  // Now that last it is set, we can call jump_to_col() to seek to the first pixel actually
  // overlapping the query. Calling jump_to_next_overlap() is required to deal with rows that are
  // not empty, but that have no pixels overlapping the query
  this->jump_to_col(_coord2->bin1_id());
  if (this->discard()) {
    this->jump_to_next_overlap();
  }
}

template <class N>
inline auto PixelSelector<N>::iterator::at_end(std::shared_ptr<const Index> index,
                                               const Dataset &pixels_bin1_id,
                                               const Dataset &pixels_bin2_id,
                                               const Dataset &pixels_count) -> iterator {
  iterator it{};
  it._index = std::move(index);
  it._bin1_id_it = pixels_bin1_id.end<std::uint64_t>();
  it._bin2_id_it = pixels_bin2_id.end<std::uint64_t>();
  it._bin2_id_last = it._bin2_id_it;
  it._count_it = pixels_count.end<N>();

  return it;
}

template <class N>
inline auto PixelSelector<N>::iterator::at_end(std::shared_ptr<const Index> index,
                                               const Dataset &pixels_bin1_id,
                                               const Dataset &pixels_bin2_id,
                                               const Dataset &pixels_count,
                                               std::shared_ptr<PixelCoordinates> coord1,
                                               std::shared_ptr<PixelCoordinates> coord2)
    -> iterator {
  if (!coord1 && !coord2) {
    return at_end(std::move(index), pixels_bin1_id, pixels_bin2_id, pixels_count);
  }
  assert(!!coord1);
  assert(!!coord2);

  assert(coord1->bin2_id() <= coord2->bin2_id());

  iterator it{};

  it._index = std::move(index);
  it._coord1 = std::move(coord1);
  it._coord2 = std::move(coord2);

  // Get the bin id for the last row overlapping the query
  auto bin1_id = it._coord1->bin2_id();

  do {
    // Get the offset to the last row overlapping the query
    auto offset = it._index->get_offset_by_bin_id(bin1_id);
    it._bin1_id_it = pixels_bin1_id.make_iterator_at_offset<std::uint64_t>(offset);
    it._bin2_id_it = pixels_bin2_id.make_iterator_at_offset<std::uint64_t>(offset);
    // Even though we do not yet know where the last iterator actually is, the jump methods expect
    // the last iterator to be set to something.
    it._bin2_id_last = pixels_bin2_id.end<std::uint64_t>();

    it._count_it = pixels_count.make_iterator_at_offset<N>(offset);

    if (offset == pixels_bin2_id.size()) {
      return it;
    }

    // Jump to the first column overlapping the query
    it.jump_to_col(it._coord2->bin1_id());

    // If discard() returns true, it means that the row corresponding to bin1_id does not contain
    // any pixel overlapping the query, so we keep looking backwards
  } while (bin1_id-- > it._coord1->bin1_id() && it.discard());

  // Now that we know the row pointed by it contains at least one pixel overlapping the query, jump
  // to the last column overlapping the query
  it.jump_to_col(it._coord2->bin2_id());

  // IMPORTANT! Do not try to set _bin2_id_last before calling it.discard(), otherwise discard will
  // always return true!
  if (it.discard()) {
    // Pixel pointed to by it does not overlap the query (i.e. matrix has 0 interactions for the
    // current row and col bin2_id()
    it._bin2_id_last = it._bin2_id_it;
  } else {
    // it points to the last pixel overlapping the query, so increment last_it and it
    it._bin2_id_last = it._bin2_id_it + 1;
    std::ignore = ++it;
  }

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
  return {PixelCoordinates{this->_index->bins_ptr(),
                           *this->_bin1_id_it,
                           *this->_bin2_id_it},
          *this->_count_it};
  // clang-format on
}

template <class N>
inline auto PixelSelector<N>::iterator::operator++() -> iterator & {
  assert(this->_bin2_id_it < this->_bin2_id_last);
  std::ignore = ++this->_bin1_id_it;
  std::ignore = ++this->_bin2_id_it;
  std::ignore = ++this->_count_it;

  if (this->discard()) {
    this->jump_to_next_overlap();
  }

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

  if (this->_bin2_id_it >= this->_bin2_id_last) {
    return;
  }

  const auto end_offset = this->_bin2_id_last.h5_offset();
  const auto row_offset = this->_index->get_offset_by_bin_id(bin_id);
  const auto current_offset = this->h5_offset();

  assert(row_offset >= current_offset);
  assert(end_offset >= current_offset);
  const auto offset = std::min(end_offset, row_offset) - current_offset;

  this->_bin1_id_it += offset;
  this->_bin2_id_it += offset;
  this->_count_it += offset;
}

template <class N>
inline void PixelSelector<N>::iterator::jump_to_col(std::uint64_t bin_id) {
  assert(this->_index);
  assert(bin_id <= this->_index->bins().size());

  if (this->_bin2_id_it >= this->_bin2_id_last) {
    return;
  }

  const auto current_row = *this->_bin1_id_it;
  const auto next_row = current_row + 1;

  const auto current_offset = conditional_static_cast<std::uint64_t>(this->h5_offset());
  const auto current_row_offset = this->_index->get_offset_by_bin_id(current_row);
  const auto next_row_offset = this->_index->get_offset_by_bin_id(next_row);
  const auto end_offset = conditional_static_cast<std::uint64_t>(this->_bin2_id_last.h5_offset());

  if (current_offset == next_row_offset) {
    return;  // Row is empty
  }

  assert(next_row_offset != 0);
  const auto row_start_offset = std::min(current_offset, current_row_offset);
  const auto row_end_offset = std::min(end_offset, next_row_offset - 1);

  if (row_start_offset == row_end_offset) {
    return;  // Row is empty
  }

  auto first = this->_bin2_id_it + (row_start_offset - current_offset);
  auto last = first + (row_end_offset - row_start_offset);
  this->_bin2_id_it = std::lower_bound(first, last, bin_id);

  const auto offset = this->_bin2_id_it.h5_offset() - current_offset;

  this->_bin1_id_it += offset;
  this->_count_it += offset;

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

template <class N>
inline void PixelSelector<N>::iterator::jump_to_next_overlap() {
  assert(this->_coord1);
  assert(this->_coord2);
  while (this->discard()) {
    // We're at/past end: return immediately
    if (this->_bin2_id_it >= this->_bin2_id_last) {
      this->jump_at_end();
      return;
    }

    const auto row = *this->_bin1_id_it;
    const auto col = *this->_bin2_id_it;
    const auto next_row = row + 1;
    const auto next_col = std::max(next_row, this->_coord2->bin1_id());

    // We may have some data left to read from the current row
    if (col < this->_coord2->bin1_id()) {
      this->jump_to_col(this->_coord2->bin1_id());
      if (!this->discard()) {
        return;
      }
    }

    // There's no more data to be read, as we're past the last column overlapping the query,
    // and the next row does not overlap the query
    if (this->_bin2_id_it >= this->_bin2_id_last || next_row > this->_coord1->bin2_id()) {
      assert(col > this->_coord2->bin2_id());
      this->jump_at_end();
      return;
    }

    this->jump(next_row, next_col);
  }
}

template <class N>
inline bool PixelSelector<N>::iterator::discard() const {
  if (!this->_coord1) {
    // Iterator is traversing the entire matrix:
    // no pixel should be discarded
    assert(!this->_coord2);
    return false;
  }

  if (this->_bin2_id_it == this->_bin2_id_last) {
    return false;
  }

  const auto overlaps_range1 = *this->_bin1_id_it >= this->_coord1->bin1_id() &&
                               *this->_bin1_id_it <= this->_coord1->bin2_id();

  const auto overlaps_range2 = *this->_bin2_id_it >= this->_coord2->bin1_id() &&
                               *this->_bin2_id_it <= this->_coord2->bin2_id();

  return !overlaps_range1 || !overlaps_range2;
}

template <class N>
inline std::size_t PixelSelector<N>::iterator::h5_offset() const noexcept {
  assert(this->_bin1_id_it.h5_offset() == this->_bin2_id_it.h5_offset());
  assert(this->_count_it.h5_offset() == this->_bin2_id_it.h5_offset());

  return this->_bin2_id_it.h5_offset();
}

template <class N>
inline void PixelSelector<N>::iterator::jump_at_end() {
  if (this->_bin2_id_it == this->_bin2_id_last) {
    return;
  }

  // Deal with iterators upstream of last
  if (this->_bin2_id_it < this->_bin2_id_last) {
    const auto offset = this->_bin2_id_last.h5_offset() - this->_bin2_id_it.h5_offset();
    this->_bin1_id_it += offset;
    this->_bin2_id_it += offset;
    this->_count_it += offset;
    return;
  }

  // Deal with iterators downstream of last
  const auto offset = this->_bin2_id_it.h5_offset() - this->_bin2_id_last.h5_offset();
  this->_bin1_id_it -= offset;
  this->_bin2_id_it -= offset;
  this->_count_it -= offset;
}

}  // namespace coolerpp
