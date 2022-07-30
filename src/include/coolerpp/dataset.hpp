// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <tsl/hopscotch_map.h>

#include <cstdint>
#include <highfive/H5DataSet.hpp>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "coolerpp/common.hpp"
#include "coolerpp/group.hpp"
#include "coolerpp/internal/generic_variant.hpp"
#include "coolerpp/internal/suppress_warnings.hpp"
#include "coolerpp/internal/variant_buff.hpp"

namespace coolerpp {

struct RootGroup;

namespace internal {
template <class T>
struct is_atomic_buffer
    : public std::disjunction<std::is_same<internal::GenericVariant, std::decay_t<T>>,
                              std::is_same<std::string, std::decay_t<T>>,
                              std::is_arithmetic<std::decay_t<T>>> {};

template <class T>
inline constexpr bool is_atomic_buffer_v = is_atomic_buffer<T>::value;
}  // namespace internal

DISABLE_WARNING_PUSH
DISABLE_WARNING_DEPRECATED_DECLARATIONS
class Dataset {
  RootGroup _root_group{};
  HighFive::DataSet _dataset{};
  mutable internal::VariantBuffer _buff{};

 public:
  template <class T>
  class iterator;
  template <class T>
  using const_iterator = iterator<T>;

  [[nodiscard]] static HighFive::DataSetCreateProps generate_default_dset_create_props(
      std::uint_fast8_t compression_lvl = DEFAULT_COMPRESSION_LEVEL,
      std::size_t chunk_size = DEFAULT_HDF5_CHUNK_SIZE);
  [[nodiscard]] static HighFive::DataSetAccessProps generate_default_dset_access_props(
      std::size_t chunk_size = DEFAULT_HDF5_CHUNK_SIZE,
      std::size_t cache_size = DEFAULT_HDF5_CACHE_SIZE);

  Dataset() = default;
  Dataset(RootGroup root_group, HighFive::DataSet dset);
  Dataset(RootGroup root_group, std::string_view path_to_dataset,
          const HighFive::DataSetAccessProps &aprops = generate_default_dset_access_props());

  template <class T, class = std::enable_if_t<std::is_arithmetic_v<T>>>
  Dataset(RootGroup root_group, std::string_view path_to_dataset, const T &type,
          std::size_t max_dim = HighFive::DataSpace::UNLIMITED,

          const HighFive::DataSetAccessProps &aprops = generate_default_dset_access_props(),
          const HighFive::DataSetCreateProps &cprops = generate_default_dset_create_props());

  Dataset(RootGroup root_group, std::string_view path_to_dataset, std::string_view longest_str,
          std::size_t max_dim = HighFive::DataSpace::UNLIMITED,
          const HighFive::DataSetAccessProps &aprops = generate_default_dset_access_props(),
          const HighFive::DataSetCreateProps &cprops = generate_default_dset_create_props());

  Dataset(const Dataset &other) = default;
  Dataset(Dataset &&other) noexcept = default;

  Dataset &operator=(const Dataset &other) = default;
  Dataset &operator=(Dataset &&other) noexcept = default;

  ~Dataset() noexcept = default;

  [[nodiscard]] std::string file_name() const;
  [[nodiscard]] std::string hdf5_path() const;
  [[nodiscard]] std::string uri() const;

  [[nodiscard]] std::size_t size() const;
  [[nodiscard]] bool empty() const;

  [[nodiscard]] HighFive::DataSet get();
  [[nodiscard]] const HighFive::DataSet &get() const;

  [[nodiscard]] RootGroup get_parent() const;

  void resize(std::size_t new_size);

  template <class N, class = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t read(std::vector<N> &buff, std::size_t num, std::size_t offset = 0) const;
  std::size_t read(std::vector<std::string> &buff, std::size_t num, std::size_t offset = 0) const;
  template <std::size_t i = 0>
  std::size_t read(internal::VariantBuffer &vbuff, std::size_t num, std::size_t offset = 0) const;

  template <class N, class = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t write(const std::vector<N> &buff, std::size_t offset = 0,
                    bool allow_dataset_resize = false);
  std::size_t write(const std::vector<std::string> &buff, std::size_t offset = 0,
                    bool allow_dataset_resize = false);
  std::size_t write(const internal::VariantBuffer &vbuff, std::size_t offset = 0,
                    bool allow_dataset_resize = false);

  template <class N, class = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t read(N &buff, std::size_t offset) const;
  std::size_t read(std::string &buff, std::size_t offset) const;
  template <std::size_t i = 0>
  std::size_t read(internal::GenericVariant &vbuff, std::size_t offset) const;

  template <class N, class = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t write(N buff, std::size_t offset = 0, bool allow_dataset_resize = false);
  std::size_t write(std::string buff, std::size_t offset = 0, bool allow_dataset_resize = false);

  std::size_t write(const internal::GenericVariant &vbuff, std::size_t offset = 0,
                    bool allow_dataset_resize = false);

  template <class BuffT, class T = remove_cvref_t<BuffT>,
            class = std::enable_if_t<!internal::is_atomic_buffer_v<T>>>
  std::size_t read_all(BuffT &buff, std::size_t offset = 0) const;

  template <class BuffT, class T = remove_cvref_t<BuffT>,
            class = std::enable_if_t<!internal::is_atomic_buffer_v<T>>>
  BuffT read_all(std::size_t offset = 0) const;

  internal::VariantBuffer read_all(std::size_t offset = 0) const;

  template <class BuffT, class T = remove_cvref_t<BuffT>,
            class = std::enable_if_t<internal::is_atomic_buffer_v<T>>>
  BuffT read(std::size_t offset) const;
  internal::GenericVariant read(std::size_t offset) const;

  template <class BuffT, class T = remove_cvref_t<BuffT>,
            class = std::enable_if_t<!internal::is_atomic_buffer_v<T>>>
  BuffT read_n(std::size_t num, std::size_t offset = 0) const;

  template <class BuffT>
  std::size_t append(const BuffT &buff);

  template <class InputIt, class UnaryOperation = identity,
            typename std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, InputIt>, int>
                * = nullptr>
  std::size_t write(InputIt first_value, InputIt last_value, std::size_t offset = 0,
                    bool allow_dataset_resize = false, UnaryOperation op = identity());

  template <class InputIt, class UnaryOperation = identity,
            typename std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, InputIt>, int>
                * = nullptr>
  std::size_t append(InputIt first_value, InputIt last_value, UnaryOperation op = identity());

  template <class BuffT>
  [[nodiscard]] BuffT read_last() const;
  [[nodiscard]] internal::GenericVariant read_last() const;

  template <class T>
  [[nodiscard]] auto begin() const -> iterator<T>;
  template <class T>
  [[nodiscard]] auto end() const -> iterator<T>;

  template <class T>
  [[nodiscard]] auto cbegin() const -> iterator<T>;
  template <class T>
  [[nodiscard]] auto cend() const -> iterator<T>;

  template <class T>
  [[nodiscard]] auto make_iterator_at_offset(std::size_t offset,
                                             std::size_t chunk_size = 64 * 1024) const
      -> iterator<T>;
  template <class T>
  [[nodiscard]] auto make_end_iterator_at_offset(std::size_t offset,
                                                 std::size_t chunk_size = 64 * 1024) const
      -> iterator<T>;

  [[nodiscard]] static std::pair<std::string, std::string> parse_uri(std::string_view uri);

 private:
  [[nodiscard]] HighFive::Selection select(std::size_t i);
  [[nodiscard]] HighFive::Selection select(std::size_t i) const;

  [[nodiscard]] HighFive::Selection select(std::size_t i1, std::size_t i2);
  [[nodiscard]] HighFive::Selection select(std::size_t i1, std::size_t i2) const;

  [[nodiscard]] static HighFive::DataSet create_fixed_str_dataset(
      RootGroup &root_grp, std::string_view path, std::size_t max_str_length, std::size_t max_dim,
      const HighFive::DataSetAccessProps &aprops, const HighFive::DataSetCreateProps &cprops);

  [[noreturn]] void throw_out_of_range_excp(std::size_t offset) const;
  [[noreturn]] void throw_out_of_range_excp(std::size_t offset, std::size_t n) const;

  [[nodiscard]] HighFive::DataType get_h5type() const;

 public:
  template <class T>
  class iterator {
    friend Dataset;

    std::shared_ptr<std::vector<T>> _buff{};
    const Dataset *_dset{};
    std::size_t _buff_capacity{};
    std::size_t _h5_chunk_start{};
    std::size_t _h5_offset{};

    explicit iterator(const Dataset &dset, std::size_t h5_offset = 0,
                      std::size_t chunk_size = 64 * 1024);

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::random_access_iterator_tag;

    iterator() = default;

    [[nodiscard]] constexpr bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] constexpr bool operator<(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator<=(const iterator &other) const noexcept;

    [[nodiscard]] constexpr bool operator>(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator>=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> value_type;
    [[nodiscard]] auto operator[](std::size_t i) const -> value_type;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
    auto operator+=(std::size_t i) -> iterator &;
    [[nodiscard]] auto operator+(std::size_t i) const -> iterator;

    auto operator--() -> iterator &;
    auto operator--(int) -> iterator;
    auto operator-=(std::size_t i) -> iterator &;
    [[nodiscard]] auto operator-(std::size_t i) const -> iterator;
    [[nodiscard]] auto operator-(const iterator &other) const -> difference_type;

   private:
    void read_chunk_at_offset(std::size_t new_offset);

    [[nodiscard]] static constexpr auto make_end_iterator(const Dataset &dset,
                                                          std::size_t chunk_size = 64 * 1024)
        -> iterator;
  };
};
DISABLE_WARNING_POP

using DatasetMap = tsl::hopscotch_map<std::string, Dataset>;

}  // namespace coolerpp

#include "../../dataset_impl.hpp"
