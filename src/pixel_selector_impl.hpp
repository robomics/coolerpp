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
                                       PixelCoordinates coords)
    : _coords(std::move(coords)),
      _index(&index),
      _pixels_bin1_id(&pixels_bin1_id),
      _pixels_bin2_id(&pixels_bin2_id),
      _pixels_count(&pixels_count),
      _empty(_coords.chrom1 == _coords.chrom2 && coords.bin1_start == coords.bin2_start) {
  if (_coords.bin2_start != 0) {
    _coords.bin2_start--;
  }
}

template <class N>
constexpr bool PixelSelector<N>::empty() const noexcept {
  return this->_empty;
}

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
  if (this->empty()) {
    return this->cend();
  }
  return iterator{*this->_index, *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                  *this->_pixels_count, this->_coords};
}

template <class N>
inline auto PixelSelector<N>::cend() const -> iterator {
  return iterator::make_end_iterator(*this->_index, *this->_pixels_bin2_id, this->_coords);
}

template <class N>
constexpr const PixelCoordinates &PixelSelector<N>::coords() const noexcept {
  return this->_coords;
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

  if (start_pos > end_pos) {
    throw std::runtime_error(fmt::format(FMT_STRING("invalid query \"{}\": query start position is "
                                                    "greater than end position ({} > {})"),
                                         query, start_pos, end_pos));
  }

  return {bins, chrom, start_pos, end_pos};
}

template <class N>
inline std::uint64_t PixelSelector<N>::compute_end_offset(const Index &index,
                                                          const Dataset &bin2_id_dset,
                                                          std::uint64_t last_bin2_id) {
  assert(last_bin2_id != 0);
  if (last_bin2_id == index.size() || last_bin2_id == index.size() + 1) {
    return bin2_id_dset.size();
  }

  const auto lower_bound = index.get_offset_by_bin_id(last_bin2_id - 1);
  const auto upper_bound = index.get_offset_by_bin_id(last_bin2_id);

  auto lb = bin2_id_dset.make_iterator_at_offset<std::uint64_t>(lower_bound, 1);
  auto ub = bin2_id_dset.make_iterator_at_offset<std::uint64_t>(upper_bound, 1);

  auto it = std::lower_bound(lb, ub, last_bin2_id);
  return lower_bound + static_cast<std::uint64_t>(std::distance(lb, it));
}

template <class N>
inline PixelSelector<N>::iterator::iterator(const Index &index, const Dataset &pixels_bin1_id,
                                            const Dataset &pixels_bin2_id,
                                            const Dataset &pixels_count,
                                            const coolerpp::PixelCoordinates &coords) noexcept
    : _index(&index),
      _first_chrom_id(index.chromosomes().get_id(coords.chrom1->name)),
      _last_chrom_id(index.chromosomes().get_id(coords.chrom2->name) + 1),
      _first_bin_id(coords.bin1_id()),
      _last_bin_id(coords.bin2_id() + 1),
      _bin1_id_it(pixels_bin1_id.make_iterator_at_offset<std::uint64_t>(
          index.get_offset_by_bin_id(_first_bin_id))),
      _bin2_id_it(pixels_bin2_id.make_iterator_at_offset<std::uint64_t>(
          index.get_offset_by_bin_id(_first_bin_id))),
      _bin2_id_last(pixels_bin2_id.make_iterator_at_offset<std::uint64_t>(
          compute_end_offset(index, pixels_bin2_id, _last_bin_id))),
      _count_it(
          pixels_count.make_iterator_at_offset<N>(index.get_offset_by_bin_id(_first_bin_id))) {
  assert(_first_chrom_id <= _last_chrom_id);
  assert(_first_bin_id < _last_bin_id);
}

template <class N>
inline auto PixelSelector<N>::iterator::make_end_iterator(const Index &index,
                                                          const Dataset &pixels_bin2_id,
                                                          const PixelCoordinates &coords)
    -> iterator {
  iterator it{};

  it._index = &index;
  it._first_chrom_id = index.chromosomes().get_id(coords.chrom1->name);
  it._last_chrom_id = index.chromosomes().get_id(coords.chrom2->name) + 1;

  it._first_bin_id = coords.bin1_id();
  it._last_bin_id = coords.bin2_id() + 1;
  it._bin2_id_last = pixels_bin2_id.make_iterator_at_offset<std::uint64_t>(
      compute_end_offset(index, pixels_bin2_id, it._last_bin_id));

  it._bin2_id_it = it._bin2_id_last;

  return it;
}

template <class N>
constexpr bool PixelSelector<N>::iterator::operator==(const iterator &other) const noexcept {
  // clang-format off
  return this->_index == other._index &&
         this->_bin2_id_it == other._bin2_id_it;
  // clang-format on
}
template <class N>
constexpr bool PixelSelector<N>::iterator::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <class N>
constexpr bool PixelSelector<N>::iterator::operator<(const iterator &other) const noexcept {
  // clang-format off
  return this->_index < other._index &&
         this->_bin2_id_it < other._bin2_id_it;
  // clang-format on
}
template <class N>
constexpr bool PixelSelector<N>::iterator::operator<=(const iterator &other) const noexcept {
  // clang-format off
  return this->_index <= other._index &&
         this->_bin2_id_it <= other._bin2_id_it;
  // clang-format on
}

template <class N>
constexpr bool PixelSelector<N>::iterator::operator>(const iterator &other) const noexcept {
  // clang-format off
  return this->_index > other._index &&
         this->_bin2_id_it > other._bin2_id_it;
  // clang-format on
}
template <class N>
constexpr bool PixelSelector<N>::iterator::operator>=(const iterator &other) const noexcept {
  // clang-format off
  return this->_index >= other._index &&
         this->_bin2_id_it >= other._bin2_id_it;
  // clang-format on
}

template <class N>
inline auto PixelSelector<N>::iterator::operator*() const -> value_type {
  // clang-format off
  assert(this->_bin2_id_it != this->_bin2_id_last);
  return {PixelCoordinates{this->_index->bins(),
                           *this->_bin1_id_it,
                           *this->_bin2_id_it},
          *this->_count_it};
  // clang-format on
}

template <class N>
inline auto PixelSelector<N>::iterator::operator++() -> iterator & {
  if (++this->_bin2_id_it >= this->_bin2_id_last) {
    ++this->_bin1_id_it;
    ++this->_count_it;
    return *this;
  }

  if (*this->_bin2_id_it == this->_last_bin_id) {
    assert(*this->_bin1_id_it < this->_last_bin_id);

    const auto new_offset = this->_index->get_offset_by_bin_id((*this->_bin1_id_it) + 1);
    const auto relative_offset = new_offset - this->_bin1_id_it.h5_offset();
    assert(relative_offset != 0);

    this->_bin1_id_it += relative_offset;
    this->_bin2_id_it += relative_offset - 1;
    this->_count_it += relative_offset;

    return *this;
  }

  std::ignore = ++this->_bin1_id_it;
  std::ignore = ++this->_count_it;

  return *this;
}

template <class N>
inline auto PixelSelector<N>::iterator::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

}  // namespace coolerpp
