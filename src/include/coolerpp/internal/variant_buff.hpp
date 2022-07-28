// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include "coolerpp/internal/generic_variant.hpp"

namespace coolerpp::internal {

class VariantBuffer {
 private:
  // Variants are listed in order from the most common to the least common for perf. reasons
  // clang-format off
  using BuffT = std::variant<
      std::vector<std::uint32_t>,
      std::vector<std::uint64_t>,
      std::vector<std::int32_t>,
      std::vector<double>,
      std::vector<std::string>,
      std::vector<std::uint8_t>,
      std::vector<std::uint16_t>,
      std::vector<std::int8_t>,
      std::vector<std::int16_t>,
      std::vector<std::int64_t>,
      std::vector<float>,
      std::vector<long double>>;
  // clang-format on
  BuffT _buff{std::vector<std::uint32_t>{}};

 public:
  VariantBuffer() = default;
  template <class T>
  explicit VariantBuffer(std::vector<T> data);
  template <class InputIt>
  VariantBuffer(InputIt first, InputIt last);
  template <class T>
  explicit VariantBuffer(std::size_t size, T default_value = 0);

  VariantBuffer(const VariantBuffer &other) = default;
  VariantBuffer(VariantBuffer &&other) noexcept = default;

  ~VariantBuffer() noexcept = default;

  VariantBuffer &operator=(const VariantBuffer &other) = default;
  VariantBuffer &operator=(VariantBuffer &&other) noexcept = default;

  template <class T>
  VariantBuffer &operator=(std::vector<T> buff) noexcept;

  template <class T>
  [[nodiscard]] typename std::vector<T>::iterator begin();
  template <class T>
  [[nodiscard]] typename std::vector<T>::iterator end();
  template <class T>
  [[nodiscard]] typename std::vector<T>::const_iterator begin() const;
  template <class T>
  typename std::vector<T>::const_iterator end() const;
  template <class T>
  [[nodiscard]] typename std::vector<T>::const_iterator cbegin() const;
  template <class T>
  [[nodiscard]] typename std::vector<T>::const_iterator cend() const;

  [[nodiscard]] std::size_t size() const noexcept;  // NOLINT(bugprone-exception-escape)
  template <class T>
  [[nodiscard]] std::size_t size() const;

  [[nodiscard]] std::size_t capacity() const noexcept;  // NOLINT(bugprone-exception-escape)
  template <class T>
  [[nodiscard]] std::size_t capacity() const;

  template <class T>
  void resize(std::size_t new_size);

  template <class T>
  void reserve(std::size_t new_size);

  [[nodiscard]] bool empty() const noexcept;  // NOLINT(bugprone-exception-escape)
  template <class T>
  [[nodiscard]] bool empty() const;

  void clear() noexcept;  // NOLINT(bugprone-exception-escape)
  template <class T>
  void clear() noexcept;

  template <class T>
  [[nodiscard]] T &at(std::size_t i);
  template <class T>
  [[nodiscard]] const T &at(std::size_t i) const;
  [[nodiscard]] GenericVariant at(std::size_t i) const;
  [[nodiscard]] GenericVariant operator[](std::size_t i) const;

  template <class T>
  [[nodiscard]] T &front();
  template <class T>
  [[nodiscard]] const T &front() const;

  template <class T>
  [[nodiscard]] T &back();
  template <class T>
  [[nodiscard]] const T &back() const;

  template <class T>
  [[nodiscard]] T *data();
  template <class T>
  [[nodiscard]] const T *data() const;

  template <class T>
  [[nodiscard]] std::vector<T> &get();
  template <class T>
  [[nodiscard]] const std::vector<T> &get() const;

  [[nodiscard]] constexpr auto get() -> BuffT &;
  [[nodiscard]] constexpr auto get() const -> const BuffT &;

  template <class T>
  [[nodiscard]] bool holds_alternative() const noexcept;
};

}  // namespace coolerpp::internal

#include "../../../variant_buff_impl.hpp"
