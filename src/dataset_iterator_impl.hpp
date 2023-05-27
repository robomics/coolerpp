// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5Exception.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Selection.hpp>
#include <memory>
#include <string>
#include <vector>

#include "coolerpp/common.hpp"
#include "coolerpp/internal/prime_number_table.hpp"
#include "coolerpp/internal/type_pretty_printer.hpp"

namespace coolerpp {

template <typename T>
inline Dataset::iterator<T>::iterator(const Dataset &dset, std::size_t h5_offset,
                                      std::size_t chunk_size, bool init)
    : _dset(&dset),
      _buff_capacity((std::min)(chunk_size, dset.size())),
      _h5_chunk_start(h5_offset),
      _h5_offset(h5_offset) {
  if (init) {
    this->read_chunk_at_offset(this->_h5_chunk_start);
  }
}

template <typename T>
constexpr bool Dataset::iterator<T>::operator==(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset == other._h5_offset;
}
template <typename T>
constexpr bool Dataset::iterator<T>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename T>
constexpr bool Dataset::iterator<T>::operator<(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset < other._h5_offset;
}
template <typename T>
constexpr bool Dataset::iterator<T>::operator<=(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset <= other._h5_offset;
}

template <typename T>
constexpr bool Dataset::iterator<T>::operator>(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset > other._h5_offset;
}
template <typename T>
constexpr bool Dataset::iterator<T>::operator>=(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset >= other._h5_offset;
}

template <typename T>
inline auto Dataset::iterator<T>::operator*() const -> value_type {
  if (!this->_buff || this->_h5_offset >= this->_h5_chunk_start + this->_buff->size()) {
    // Read first chunk
    this->read_chunk_at_offset(this->_h5_offset);
  } else if (this->_h5_offset - this->_h5_chunk_start >= this->_buff->size()) {
    // Iterator was decremented one or more times since the last dereference, thus we assume the
    // iterator is being used to traverse the dataset backward
    this->_h5_chunk_start =
        this->_h5_offset - (std::min)(this->_buff->size() - 1, this->_h5_offset);
    this->read_chunk_at_offset(this->_h5_chunk_start);
  }

  assert(this->_buff);
  assert(this->_dset);
  assert(this->_h5_offset < this->_dset->size());
  assert(this->_h5_chunk_start <= this->_h5_offset);
  assert(this->_h5_offset - this->_h5_chunk_start < this->_buff->size());
  return (*this->_buff)[this->_h5_offset - this->_h5_chunk_start];
}

template <typename T>
inline auto Dataset::iterator<T>::operator[](std::size_t i) const -> value_type {
  return *(*this + i);
}

template <typename T>
inline auto Dataset::iterator<T>::operator++() -> iterator & {
  return (*this) += 1;
}

template <typename T>
inline auto Dataset::iterator<T>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  if (this->_h5_offset > this->_h5_chunk_start + this->_buff_capacity) {
    this->read_chunk_at_offset(this->_h5_offset);
  }
  return it;
}

template <typename T>
inline auto Dataset::iterator<T>::operator+=(std::size_t i) -> iterator & {
  assert(this->_dset);
  assert(this->_h5_offset + i <= this->_dset->size());
  this->_h5_offset += i;
  return *this;
}

template <typename T>
inline auto Dataset::iterator<T>::operator+(std::size_t i) const -> iterator {
  assert(this->_dset);
  assert(this->_buff);
  const auto new_offset = this->_h5_offset + i;
  assert(new_offset <= this->_dset->size());

  if (!this->_buff || this->_h5_chunk_start + this->_buff->size() < new_offset) {
    return iterator(*this->_dset, new_offset);
  }

  auto it = *this;
  return it += i;
}

template <typename T>
inline auto Dataset::iterator<T>::operator--() -> iterator & {
  assert(this->_h5_offset != 0);
  return (*this) -= 1;
}

template <typename T>
inline auto Dataset::iterator<T>::operator--(int) -> iterator {
  auto it = *this;
  std::ignore = --(*this);
  if (this->_h5_offset < this->_h5_chunk_start) {
    this->read_chunk_at_offset(this->_h5_offset -
                               (std::min)(this->_buff_capacity - 1, this->_h5_offset));
  }
  return it;
}

template <typename T>
inline auto Dataset::iterator<T>::operator-=(std::size_t i) -> iterator & {
  assert(this->_h5_offset >= i);
  this->_h5_offset -= i;
  return *this;
}

template <typename T>
inline auto Dataset::iterator<T>::operator-(std::size_t i) const -> iterator {
  assert(this->_h5_offset >= i);
  const auto new_offset = this->_h5_offset - i;
  if (new_offset >= this->_h5_chunk_start) {
    auto it = *this;
    return it -= i;
  }

  assert(this->_dset);
  return iterator(*this->_dset, new_offset);
}

template <typename T>
inline auto Dataset::iterator<T>::operator-(const iterator &other) const -> difference_type {
  return static_cast<difference_type>(this->_h5_offset) -
         static_cast<difference_type>(other._h5_offset);
}

template <typename T>
constexpr std::uint64_t Dataset::iterator<T>::h5_offset() const noexcept {
  return this->_h5_offset;
}

template <typename T>
constexpr std::size_t Dataset::iterator<T>::underlying_buff_capacity() const noexcept {
  return this->_buff_capacity;
}

template <typename T>
constexpr const Dataset &Dataset::iterator<T>::dataset() const noexcept {
  return *this->_dset;
}

template <typename T>
inline void Dataset::iterator<T>::read_chunk_at_offset(std::size_t new_offset) const {
  assert(this->_dset);

  if (new_offset == this->_dset->size()) {
    this->_buff = nullptr;
    this->_h5_chunk_start = this->_dset->size();
    return;
  }

  if (!this->_buff || !this->_buff.unique()) {
    //  This should be fine, as copying Dataset::iterator is not thread-safe anyway
    this->_buff = std::make_shared<std::vector<T>>(this->_buff_capacity);
  }

  const auto buff_size = (std::min)(this->_buff_capacity, this->_dset->size() - new_offset);
  this->_buff->resize(buff_size);
  this->_dset->read(*this->_buff, buff_size, new_offset);

  this->_h5_chunk_start = new_offset;
}

template <typename T>
constexpr auto Dataset::iterator<T>::make_end_iterator(const Dataset &dset, std::size_t chunk_size)
    -> iterator {
  iterator it{};
  it._buff = nullptr;
  it._buff_capacity = (std::min)(chunk_size, dset.size());
  it._dset = &dset;
  it._h5_offset = dset.size();
  it._h5_chunk_start = it._h5_offset;

  return it;
}

}  // namespace coolerpp
