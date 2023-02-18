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
#include <vector>

#include "coolerpp/common.hpp"
#include "coolerpp/internal/type_pretty_printer.hpp"

namespace coolerpp {

template <class T, class>
inline Dataset::Dataset(RootGroup root_group, std::string_view path_to_dataset,
                        [[maybe_unused]] const T &type, std::size_t max_dim,
                        const HighFive::DataSetAccessProps &aprops,
                        const HighFive::DataSetCreateProps &cprops)
    : Dataset(root_group,
              root_group().createDataSet<T>(std::string{path_to_dataset},
                                            HighFive::DataSpace({0}, {max_dim}), cprops, aprops)) {}

template <class N, class>
inline std::size_t Dataset::read(std::vector<N> &buff, std::size_t num, std::size_t offset) const {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  if (offset + num > this->size()) {
    this->throw_out_of_range_excp(offset, num);
  }

  auto h5type = this->get_h5type();
  buff.resize(num);
  this->select(offset, num).read(buff.data(), HighFive::create_datatype<N>());

  return offset + num;
}

template <std::size_t i>
inline std::size_t Dataset::read(internal::GenericVariant &vbuff, std::size_t offset) const {
  if constexpr (i == 0) {
    if (offset >= this->size()) {
      this->throw_out_of_range_excp(offset);
    }
  }

  using VBuffT = internal::GenericVariant;
  if constexpr (i < std::variant_size_v<VBuffT>) {
    using T = std::variant_alternative_t<i, VBuffT>;

    const auto h5type = this->get_h5type();
    if constexpr (is_string_v<T>) {
      if (h5type.isFixedLenStr() || h5type.isVariableStr()) {
        goto READ_VARIANT;  // NOLINT
      }
    }
    if (h5type != HighFive::create_datatype<T>()) {
      return read<i + 1>(vbuff, offset);
    }

#if defined(__GNUC__) && !defined(__clang__)
    // Workaround for buggy -Wunused-label on GCC
    goto READ_VARIANT;  // NOLINT
#endif

  READ_VARIANT:
    if (!std::holds_alternative<T>(vbuff)) {
      vbuff = T{};
    }

    return this->read(std::get<T>(vbuff), offset);
  }

  unreachable_code();
}

template <class BuffT, class T, class>
inline std::size_t Dataset::read_all(BuffT &buff, std::size_t offset) const {
  const auto num = offset > this->size() ? std::size_t(0) : this->size() - offset;
  return this->read(buff, num, offset);
}

template <class BuffT, class T, class>
inline BuffT Dataset::read_all(std::size_t offset) const {
  BuffT buff{};
  this->read_all(buff, offset);
  return buff;
}

template <class N, class>
inline std::size_t Dataset::write(const std::vector<N> &buff, std::size_t offset,
                                  bool allow_dataset_resize) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  if (offset + buff.size() > this->size()) {
    if (allow_dataset_resize) {
      this->resize(offset + buff.size());
    } else {
      this->throw_out_of_range_excp(offset, buff.size());
    }
  }

  this->select(offset, buff.size()).write(buff);
  return offset + buff.size();
}

template <class N, class>
inline std::size_t Dataset::read(N &buff, std::size_t offset) const {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  if (offset >= this->size()) {
    this->throw_out_of_range_excp(offset);
  }

  this->select(offset, 1).read(&buff, HighFive::create_datatype<N>());
  return offset + 1;
}

template <std::size_t i>
inline std::size_t Dataset::read(internal::VariantBuffer &vbuff, std::size_t num,
                                 std::size_t offset) const {
  if constexpr (i == 0) {
    if (offset + num > this->size()) {
      this->throw_out_of_range_excp(offset, num);
    }
  }

  using VBuffT = remove_cvref_t<decltype(vbuff.get())>;
  if constexpr (i < std::variant_size_v<VBuffT>) {
    using VT = std::variant_alternative_t<i, VBuffT>;
    using T = typename VT::value_type;

    auto h5type = this->get_h5type();
    if constexpr (is_string_v<T>) {
      if (h5type.isFixedLenStr() || h5type.isVariableStr()) {
        goto READ_VARIANT;  // NOLINT
      }
    }
    if (h5type != HighFive::create_datatype<T>()) {
      return read<i + 1>(vbuff, num, offset);
    }

#if defined(__GNUC__) && !defined(__clang__)
    // Workaround for buggy -Wunused-label on GCC
    goto READ_VARIANT;  // NOLINT
#endif

  READ_VARIANT:
    if (!std::holds_alternative<VT>(vbuff.get())) {
      vbuff = std::vector<T>(num);
    }

    return this->read(vbuff.get<T>(), num, offset);
  }

  unreachable_code();
}

template <class N, class>
inline std::size_t Dataset::write(N buff, std::size_t offset, bool allow_dataset_resize) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  if (offset >= this->size()) {
    if (allow_dataset_resize) {
      this->resize(offset + 1);
    } else {
      this->throw_out_of_range_excp(offset);
    }
  }

  this->select(offset).write(buff);
  return offset + 1;
}

template <class BuffT, class T, class>
inline BuffT Dataset::read(std::size_t offset) const {
  BuffT buff{};
  this->read(buff, offset);
  return buff;
}

template <class BuffT, class T, class>
inline BuffT Dataset::read_n(std::size_t num, std::size_t offset) const {
  BuffT buff{num};
  this->read(buff, offset);
  return buff;
}

template <class BuffT>
inline std::size_t Dataset::append(const BuffT &buff) {
  return this->write(buff, this->size(), true);
}

template <class InputIt, class UnaryOperation,
          typename std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, InputIt>, int> *>
inline std::size_t Dataset::write(InputIt first_value, InputIt last_value, std::size_t offset,
                                  bool allow_dataset_resize, UnaryOperation op) {
  using T = remove_cvref_t<decltype(op(*first_value))>;
  constexpr std::size_t buffer_capacity = is_string_v<T> ? 256 : (64 * 1024 * 1024) / sizeof(T);
  if (this->_buff.holds_alternative<T>()) {
    this->_buff.resize<T>(buffer_capacity);
  } else {
    this->_buff = std::vector<T>(buffer_capacity);
  }

  auto &buff = this->_buff.get<T>();
  buff.clear();

  while (first_value != last_value) {
    if (buff.size() == buff.capacity()) {
      this->write(buff, offset, allow_dataset_resize);
      offset += buff.size();
      buff.clear();
    }

    buff.emplace_back(op(*first_value++));
  }

  if (!buff.empty()) {
    this->write(buff, offset, allow_dataset_resize);
    offset += buff.size();
  }

  return offset;
}

template <class InputIt, class UnaryOperation,
          typename std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, InputIt>, int> *>
inline std::size_t Dataset::append(InputIt first_value, InputIt last_value, UnaryOperation op) {
  return this->write(first_value, last_value, this->size(), true, op);
}

template <class BuffT>
inline BuffT Dataset::read_last() const {
  if (this->empty()) {
    this->throw_out_of_range_excp(0);
  }
  BuffT buff{};
  this->read(buff, this->size() - 1);

  return buff;
}

template <class T>
inline void Dataset::write_attribute(std::string_view key, const T &value,
                                     bool overwrite_if_exists) {
  Attribute::write(this->_dataset, key, value, overwrite_if_exists);
}

template <class T>
inline T Dataset::read_attribute(std::string_view key) const {
  return Attribute::read<T>(this->_dataset, key);
}

template <class T>
inline auto Dataset::begin() const -> iterator<T> {
  return iterator<T>(*this);
}

template <class T>
inline auto Dataset::end() const -> iterator<T> {
  return iterator<T>::make_end_iterator(*this);
}

template <class T>
inline auto Dataset::cbegin() const -> iterator<T> {
  return this->begin<T>();
}

template <class T>
inline auto Dataset::cend() const -> iterator<T> {
  return this->end<T>();
}

template <class T>
inline auto Dataset::make_iterator_at_offset(std::size_t offset, std::size_t chunk_size) const
    -> iterator<T> {
  return iterator<T>(*this, offset, chunk_size);
}

template <class T>
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

template <class T>
constexpr bool Dataset::iterator<T>::operator==(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset == other._h5_offset;
}
template <class T>
constexpr bool Dataset::iterator<T>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <class T>
constexpr bool Dataset::iterator<T>::operator<(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset < other._h5_offset;
}
template <class T>
constexpr bool Dataset::iterator<T>::operator<=(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset <= other._h5_offset;
}

template <class T>
constexpr bool Dataset::iterator<T>::operator>(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset > other._h5_offset;
}
template <class T>
constexpr bool Dataset::iterator<T>::operator>=(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset >= other._h5_offset;
}

template <class T>
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

template <class T>
inline auto Dataset::iterator<T>::operator[](std::size_t i) const -> value_type {
  return *(*this + i);
}

template <class T>
inline auto Dataset::iterator<T>::operator++() -> iterator & {
  return (*this) += 1;
}

template <class T>
inline auto Dataset::iterator<T>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  if (this->_h5_offset > this->_h5_chunk_start + this->_buff_capacity) {
    this->read_chunk_at_offset(this->_h5_offset);
  }
  return it;
}

template <class T>
inline auto Dataset::iterator<T>::operator+=(std::size_t i) -> iterator & {
  assert(this->_dset);
  assert(this->_h5_offset + i <= this->_dset->size());
  this->_h5_offset += i;
  return *this;
}

template <class T>
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

template <class T>
inline auto Dataset::iterator<T>::operator--() -> iterator & {
  assert(this->_h5_offset != 0);
  return (*this) -= 1;
}

template <class T>
inline auto Dataset::iterator<T>::operator--(int) -> iterator {
  auto it = *this;
  std::ignore = --(*this);
  if (this->_h5_offset < this->_h5_chunk_start) {
    this->read_chunk_at_offset(this->_h5_offset -
                               (std::min)(this->_buff_capacity - 1, this->_h5_offset));
  }
  return it;
}

template <class T>
inline auto Dataset::iterator<T>::operator-=(std::size_t i) -> iterator & {
  assert(this->_h5_offset >= i);
  this->_h5_offset -= i;
  return *this;
}

template <class T>
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

template <class T>
inline auto Dataset::iterator<T>::operator-(const iterator &other) const -> difference_type {
  return static_cast<difference_type>(this->_h5_offset) -
         static_cast<difference_type>(other._h5_offset);
}

template <class T>
constexpr std::uint64_t Dataset::iterator<T>::h5_offset() const noexcept {
  return this->_h5_offset;
}

template <class T>
constexpr std::size_t Dataset::iterator<T>::underlying_buff_capacity() const noexcept {
  return this->_buff_capacity;
}

template <class T>
constexpr const Dataset &Dataset::iterator<T>::dataset() const noexcept {
  return *this->_dset;
}

template <class T>
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

template <class T>
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
