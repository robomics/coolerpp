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

template <typename N, std::size_t CHUNK_SIZE>
inline PixelSelector<N, CHUNK_SIZE>::PixelSelector(std::shared_ptr<const Index> index,
                                                   const Dataset &pixels_bin1_id,
                                                   const Dataset &pixels_bin2_id,
                                                   const Dataset &pixels_count,
                                                   PixelCoordinates coords) noexcept
    : PixelSelector(std::move(index), pixels_bin1_id, pixels_bin2_id, pixels_count, coords,
                    std::move(coords)) {}

template <typename N, std::size_t CHUNK_SIZE>
inline PixelSelector<N, CHUNK_SIZE>::PixelSelector(std::shared_ptr<const Index> index,
                                                   const Dataset &pixels_bin1_id,
                                                   const Dataset &pixels_bin2_id,
                                                   const Dataset &pixels_count,
                                                   PixelCoordinates coord1,
                                                   PixelCoordinates coord2) noexcept
    : _coord1(std::move(coord1)),
      _coord2(std::move(coord2)),
      _index(std::move(index)),
      _pixels_bin1_id(&pixels_bin1_id),
      _pixels_bin2_id(&pixels_bin2_id),
      _pixels_count(&pixels_count) {
  assert(_index);
}

template <typename N, std::size_t CHUNK_SIZE>
inline PixelSelector<N, CHUNK_SIZE>::PixelSelector(std::shared_ptr<const Index> index,
                                                   const Dataset &pixels_bin1_id,
                                                   const Dataset &pixels_bin2_id,
                                                   const Dataset &pixels_count) noexcept
    : PixelSelector<N, CHUNK_SIZE>(std::move(index), pixels_bin1_id, pixels_bin2_id, pixels_count,
                                   PixelCoordinates{}, PixelCoordinates{}) {}

template <typename N, std::size_t CHUNK_SIZE>
template <std::size_t CHUNK_SIZE_OTHER>
inline bool PixelSelector<N, CHUNK_SIZE>::operator==(
    const PixelSelector<N, CHUNK_SIZE_OTHER> &other) const noexcept {
  return this->begin() == other.begin() && this->end() == other.end();
}

template <typename N, std::size_t CHUNK_SIZE>
template <std::size_t CHUNK_SIZE_OTHER>
inline bool PixelSelector<N, CHUNK_SIZE>::operator!=(
    const PixelSelector<N, CHUNK_SIZE_OTHER> &other) const noexcept {
  return !(*this == other);
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto PixelSelector<N, CHUNK_SIZE>::begin() const -> iterator {
  return this->cbegin();
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto PixelSelector<N, CHUNK_SIZE>::end() const -> iterator {
  return this->cend();
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto PixelSelector<N, CHUNK_SIZE>::cbegin() const -> iterator {
  if (!this->_coord1) {
    assert(!this->_coord2);
    return iterator{this->_index, *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                    *this->_pixels_count};
  }

  return iterator{this->_index,         *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                  *this->_pixels_count, this->_coord1,          this->_coord2};
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto PixelSelector<N, CHUNK_SIZE>::cend() const -> iterator {
  if (!this->_coord1) {
    assert(!this->_coord2);
    return iterator::at_end(this->_index, *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                            *this->_pixels_count);
  }
  return iterator::at_end(this->_index, *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                          *this->_pixels_count, this->_coord1, this->_coord2);
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr const PixelCoordinates &PixelSelector<N, CHUNK_SIZE>::coord1() const noexcept {
  assert(this->_coord1);
  return *this->_coord1;
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr const PixelCoordinates &PixelSelector<N, CHUNK_SIZE>::coord2() const noexcept {
  assert(this->_coord1);
  return *this->_coord2;
}

template <typename N, std::size_t CHUNK_SIZE>
inline PixelSelector<N, CHUNK_SIZE>::iterator::iterator(std::shared_ptr<const Index> index,
                                                        const Dataset &pixels_bin1_id,
                                                        const Dataset &pixels_bin2_id,
                                                        const Dataset &pixels_count)
    : _bin1_id_it(pixels_bin1_id.begin<std::uint64_t, CHUNK_SIZE>()),
      _bin2_id_it(pixels_bin2_id.begin<std::uint64_t, CHUNK_SIZE>()),
      _count_it(pixels_count.begin<N, CHUNK_SIZE>()),
      _index(std::move(index)),
      _coord1(nullptr),
      _coord2(nullptr),
      _h5_end_offset(pixels_bin2_id.size()) {}

template <typename N, std::size_t CHUNK_SIZE>
inline PixelSelector<N, CHUNK_SIZE>::iterator::iterator(std::shared_ptr<const Index> index,
                                                        const Dataset &pixels_bin1_id,
                                                        const Dataset &pixels_bin2_id,
                                                        const Dataset &pixels_count,
                                                        PixelCoordinates coord1,
                                                        PixelCoordinates coord2)
    : _index(std::move(index)),
      _coord1(std::make_shared<PixelCoordinates>(std::move(coord1))),
      _coord2(std::make_shared<PixelCoordinates>(std::move(coord2))) {
  assert(_coord1);
  assert(_coord2);
  assert(_coord1->bin1.id() <= _coord1->bin2.id());
  assert(_coord2->bin1.id() <= _coord2->bin2.id());

  // Set iterator to the first row overlapping the query (i.e. the first bin overlapping coord1)
  auto offset = _index->get_offset_by_bin_id(_coord1->bin1.id());
  _bin1_id_it = pixels_bin1_id.make_iterator_at_offset<std::uint64_t, CHUNK_SIZE>(offset);
  _bin2_id_it = pixels_bin2_id.make_iterator_at_offset<std::uint64_t, CHUNK_SIZE>(offset);
  _count_it = pixels_count.make_iterator_at_offset<N, CHUNK_SIZE>(offset);

  // Init last iterator. at_end will return an iterator pointing to the pixel past the last pixel
  // overlapping the query
  auto it = iterator::at_end(this->_index, pixels_bin1_id, pixels_bin2_id, pixels_count, _coord1,
                             _coord2);
  _h5_end_offset = it._h5_end_offset;
  assert(!this->is_past_end());

  if (_bin2_id_it.h5_offset() == pixels_bin2_id.size()) {
    return;
  }

  // Now that last it is set, we can call jump_to_col() to seek to the first pixel actually
  // overlapping the query. Calling jump_to_next_overlap() is required to deal with rows that are
  // not empty, but that have no pixels overlapping the query
  this->jump_to_col(_coord2->bin1.id());
  if (this->discard()) {
    this->jump_to_next_overlap();
  }
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto PixelSelector<N, CHUNK_SIZE>::iterator::at_end(std::shared_ptr<const Index> index,
                                                           const Dataset &pixels_bin1_id,
                                                           const Dataset &pixels_bin2_id,
                                                           const Dataset &pixels_count)
    -> iterator {
  iterator it{};
  it._index = std::move(index);
  it._bin1_id_it = pixels_bin1_id.end<std::uint64_t, CHUNK_SIZE>();
  it._bin2_id_it = pixels_bin2_id.end<std::uint64_t, CHUNK_SIZE>();
  it._h5_end_offset = it._bin2_id_it.h5_offset();
  it._count_it = pixels_count.end<N, CHUNK_SIZE>();

  return it;
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto PixelSelector<N, CHUNK_SIZE>::iterator::at_end(std::shared_ptr<const Index> index,
                                                           const Dataset &pixels_bin1_id,
                                                           const Dataset &pixels_bin2_id,
                                                           const Dataset &pixels_count,
                                                           PixelCoordinates coord1,
                                                           PixelCoordinates coord2) -> iterator {
  std::shared_ptr<const PixelCoordinates> coord1_ =
      !!coord1 ? std::make_shared<const PixelCoordinates>(std::move(coord1)) : nullptr;
  std::shared_ptr<const PixelCoordinates> coord2_ =
      !!coord2 ? std::make_shared<const PixelCoordinates>(std::move(coord2)) : nullptr;

  return PixelSelector<N, CHUNK_SIZE>::iterator::at_end(std::move(index), pixels_bin1_id,
                                                        pixels_bin2_id, pixels_count,
                                                        std::move(coord1_), std::move(coord2_));
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto PixelSelector<N, CHUNK_SIZE>::iterator::at_end(
    std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
    const Dataset &pixels_bin2_id, const Dataset &pixels_count,
    std::shared_ptr<const PixelCoordinates> coord1, std::shared_ptr<const PixelCoordinates> coord2)
    -> iterator {
  if (!coord1 && !coord2) {
    return at_end(std::move(index), pixels_bin1_id, pixels_bin2_id, pixels_count);
  }
  assert(!!coord1);
  assert(!!coord2);

  assert(coord1->bin2.id() <= coord2->bin2.id());

  iterator it{};

  it._index = std::move(index);
  it._coord1 = std::move(coord1);
  it._coord2 = std::move(coord2);

  // Get the bin id for the last row overlapping the query
  auto bin1_id = it._coord1->bin2.id();

  if (bin1_id != 0) {
    // Efficiently deal with the possiblity that a cooler file does not have any
    // interactions overlapping or upstream of the queried region
    const auto offset = it._index->get_offset_by_bin_id(bin1_id);
    if (offset == 0) {
      it._bin1_id_it = pixels_bin1_id.begin<std::uint64_t, CHUNK_SIZE>();
      it._bin2_id_it = pixels_bin2_id.begin<std::uint64_t, CHUNK_SIZE>();
      it._count_it = pixels_count.begin<N, CHUNK_SIZE>();
      it._h5_end_offset = it._bin2_id_it.h5_offset();
      return it;
    }
  }

  do {
    // Get the offset to the last row overlapping the query
    const auto offset = it._index->get_offset_by_bin_id(bin1_id);
    it._bin1_id_it = pixels_bin1_id.make_iterator_at_offset<std::uint64_t, CHUNK_SIZE>(offset);
    it._bin2_id_it = pixels_bin2_id.make_iterator_at_offset<std::uint64_t, CHUNK_SIZE>(offset);
    // Even though we do not yet know where the last iterator actually is, the jump methods expect
    // _h5_end_offset to be set to something.
    it._h5_end_offset = pixels_bin2_id.size();

    it._count_it = pixels_count.make_iterator_at_offset<N, CHUNK_SIZE>(offset);

    if (offset == pixels_bin2_id.size()) {
      return it;
    }

    // Jump to the first column overlapping the query
    it.jump_to_col(it._coord2->bin1.id());

    // If discard() returns true, it means that the row corresponding to bin1_id does not contain
    // any pixel overlapping the query, so we keep looking backwards
  } while (bin1_id-- > it._coord1->bin1.id() && it.discard());

  // Now that we know the row pointed by it contains at least one pixel overlapping the query, jump
  // to the last column overlapping the query
  it.jump_to_col(it._coord2->bin2.id());

  // IMPORTANT! Do not try to set _bin2_id_last before calling it.discard(), otherwise discard will
  // always return true!
  if (it.discard()) {
    // Pixel pointed to by it does not overlap the query (i.e. matrix has 0 interactions for the
    // current row and col bin2_id()
    it._h5_end_offset = it._bin2_id_it.h5_offset();
  } else {
    // it points to the last pixel overlapping the query, so increment last_it and it
    it._h5_end_offset = it._bin2_id_it.h5_offset() + 1;
    std::ignore = ++it;
  }

  return it;
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool PixelSelector<N, CHUNK_SIZE>::iterator::operator==(
    const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it == other._bin2_id_it;
}
template <typename N, std::size_t CHUNK_SIZE>
constexpr bool PixelSelector<N, CHUNK_SIZE>::iterator::operator!=(
    const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool PixelSelector<N, CHUNK_SIZE>::iterator::operator<(
    const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it < other._bin2_id_it;
}
template <typename N, std::size_t CHUNK_SIZE>
constexpr bool PixelSelector<N, CHUNK_SIZE>::iterator::operator<=(
    const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it <= other._bin2_id_it;
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool PixelSelector<N, CHUNK_SIZE>::iterator::operator>(
    const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it > other._bin2_id_it;
}
template <typename N, std::size_t CHUNK_SIZE>
constexpr bool PixelSelector<N, CHUNK_SIZE>::iterator::operator>=(
    const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it >= other._bin2_id_it;
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto PixelSelector<N, CHUNK_SIZE>::iterator::operator*() const -> const_reference {
  assert(this->_index);
  assert_within_bound();
  if (this->pixel_is_outdated()) {
    this->read_pixel();
  }
  return this->_value;
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto PixelSelector<N, CHUNK_SIZE>::iterator::operator->() const -> const_pointer {
  return &(*(*this));
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto PixelSelector<N, CHUNK_SIZE>::iterator::operator++() -> iterator & {
  assert_within_bound();
  std::ignore = ++this->_bin1_id_it;
  std::ignore = ++this->_bin2_id_it;
  std::ignore = ++this->_count_it;

  if (this->discard()) {
    this->jump_to_next_overlap();
  }

  // signal _value is outdated
  this->_value.count = 0;

  return *this;
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto PixelSelector<N, CHUNK_SIZE>::iterator::operator++(int) -> iterator {
  if (this->_bin1_id_it.underlying_buff_num_available_fwd() <= 1) {
    this->refresh();
  }
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <typename N, std::size_t CHUNK_SIZE>
inline void PixelSelector<N, CHUNK_SIZE>::iterator::jump_to_row(std::uint64_t bin_id) {
  assert(this->_index);
  assert(bin_id <= this->_index->bins().size());

  if (this->is_at_end() || this->is_past_end()) {
    return;
  }

  const auto end_offset = this->_h5_end_offset;
  const auto row_offset = this->_index->get_offset_by_bin_id(bin_id);
  const auto current_offset = this->h5_offset();

  assert(row_offset >= current_offset);
  assert(end_offset >= current_offset);
  const auto offset = (std::min)(end_offset, row_offset) - current_offset;

  this->_bin1_id_it += offset;
  this->_bin2_id_it += offset;
  this->_count_it += offset;
}

template <typename N, std::size_t CHUNK_SIZE>
inline void PixelSelector<N, CHUNK_SIZE>::iterator::jump_to_col(std::uint64_t bin_id) {
  assert(this->_index);
  assert(bin_id <= this->_index->bins().size());

  if (this->is_at_end() || this->is_past_end()) {
    return;
  }

  const auto current_row = *this->_bin1_id_it;
  const auto next_row = current_row + 1;

  const auto current_offset = conditional_static_cast<std::uint64_t>(this->h5_offset());
  const auto current_row_offset = this->_index->get_offset_by_bin_id(current_row);
  const auto next_row_offset = this->_index->get_offset_by_bin_id(next_row);
  const auto end_offset = this->_h5_end_offset;

  if (current_offset == next_row_offset) {
    return;  // Row is empty
  }

  assert(next_row_offset != 0);
  const auto row_start_offset = (std::min)(current_offset, current_row_offset);
  const auto row_end_offset = (std::min)(end_offset, next_row_offset - 1);

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
  assert(!this->is_past_end());
}

template <typename N, std::size_t CHUNK_SIZE>
inline void PixelSelector<N, CHUNK_SIZE>::iterator::jump(std::uint64_t bin1_id,
                                                         std::uint64_t bin2_id) {
  assert(bin1_id <= bin2_id);

  this->jump_to_row(bin1_id);
  if (bin2_id != bin1_id) {
    this->jump_to_col(bin2_id);
  }
}

template <typename N, std::size_t CHUNK_SIZE>
inline void PixelSelector<N, CHUNK_SIZE>::iterator::jump_to_next_overlap() {
  assert(this->_coord1);
  assert(this->_coord2);
  while (this->discard()) {
    // We're at/past end: return immediately
    if (this->is_at_end() || this->is_past_end()) {
      this->jump_at_end();
      return;
    }

    const auto row = *this->_bin1_id_it;
    const auto col = *this->_bin2_id_it;
    const auto next_row = row + 1;
    const auto next_col = (std::max)(next_row, this->_coord2->bin1.id());

    // We may have some data left to read from the current row
    if (col < this->_coord2->bin1.id()) {
      this->jump_to_col(this->_coord2->bin1.id());
      if (!this->discard()) {
        return;
      }
    }

    // There's no more data to be read, as we're past the last column overlapping the query,
    // and the next row does not overlap the query
    if (this->is_at_end() || this->is_past_end() || next_row > this->_coord1->bin2.id()) {
      // assert(col > this->_coord2->bin2.id());  // This is not always true for trans queries
      this->jump_at_end();
      return;
    }

    this->jump(next_row, next_col);
  }
}

template <typename N, std::size_t CHUNK_SIZE>
inline bool PixelSelector<N, CHUNK_SIZE>::iterator::discard() const {
  if (!this->_coord1) {
    // Iterator is traversing the entire matrix:
    // no pixel should be discarded
    assert(!this->_coord2);
    return false;
  }

  if (this->is_at_end()) {
    return false;
  }

  const auto overlaps_range1 = *this->_bin1_id_it >= this->_coord1->bin1.id() &&
                               *this->_bin1_id_it <= this->_coord1->bin2.id();

  const auto overlaps_range2 = *this->_bin2_id_it >= this->_coord2->bin1.id() &&
                               *this->_bin2_id_it <= this->_coord2->bin2.id();

  return !overlaps_range1 || !overlaps_range2;
}

template <typename N, std::size_t CHUNK_SIZE>
inline std::size_t PixelSelector<N, CHUNK_SIZE>::iterator::h5_offset() const noexcept {
  assert(this->_bin1_id_it.h5_offset() == this->_bin2_id_it.h5_offset());
  assert(this->_count_it.h5_offset() == this->_bin2_id_it.h5_offset());

  return this->_bin2_id_it.h5_offset();
}

template <typename N, std::size_t CHUNK_SIZE>
inline void PixelSelector<N, CHUNK_SIZE>::iterator::jump_at_end() {
  if (this->is_at_end()) {
    return;
  }

  // Deal with iterators upstream of last
  if (!this->is_past_end()) {
    const auto offset = this->_h5_end_offset - this->_bin2_id_it.h5_offset();
    this->_bin1_id_it += offset;
    this->_bin2_id_it += offset;
    this->_count_it += offset;
    return;
  }

  // Deal with iterators downstream of last
  const auto offset = this->_bin2_id_it.h5_offset() - this->_h5_end_offset;
  this->_bin1_id_it -= offset;
  this->_bin2_id_it -= offset;
  this->_count_it -= offset;
}

template <typename N, std::size_t CHUNK_SIZE>
inline void PixelSelector<N, CHUNK_SIZE>::iterator::refresh() {
  const auto h5_offset = this->_bin1_id_it.h5_offset();

  const auto &bin1_dset = this->_bin1_id_it.dataset();
  const auto &bin2_dset = this->_bin2_id_it.dataset();
  const auto &count_dset = this->_count_it.dataset();

  this->_bin1_id_it =
      bin1_dset.template make_iterator_at_offset<std::uint64_t, CHUNK_SIZE>(h5_offset);
  this->_bin2_id_it =
      bin2_dset.template make_iterator_at_offset<std::uint64_t, CHUNK_SIZE>(h5_offset);
  this->_count_it = count_dset.template make_iterator_at_offset<N, CHUNK_SIZE>(h5_offset);
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool PixelSelector<N, CHUNK_SIZE>::iterator::pixel_is_outdated() const noexcept {
  return this->_value.count == 0;
}

template <typename N, std::size_t CHUNK_SIZE>
inline void PixelSelector<N, CHUNK_SIZE>::iterator::read_pixel() const {
  assert(this->pixel_is_outdated());
  this->_value.coords = PixelCoordinates{this->_index->bins().at(*this->_bin1_id_it),
                                         this->_index->bins().at(*this->_bin2_id_it)};
  this->_value.count = *this->_count_it;
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr void PixelSelector<N, CHUNK_SIZE>::iterator::assert_within_bound() const noexcept {
  assert(this->_bin2_id_it.h5_offset() < this->_h5_end_offset);
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool PixelSelector<N, CHUNK_SIZE>::iterator::is_at_end() const noexcept {
  return this->_bin2_id_it.h5_offset() == this->_h5_end_offset;
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool PixelSelector<N, CHUNK_SIZE>::iterator::is_past_end() const noexcept {
  return this->_bin2_id_it.h5_offset() > this->_h5_end_offset;
}
}  // namespace coolerpp
