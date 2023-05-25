// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>

#include "coolerpp/bin_table.hpp"
#include "coolerpp/dataset.hpp"
#include "coolerpp/pixel_selector.hpp"

namespace coolerpp {

class Weights {
 public:
  enum class Type { INFER, DIVISIVE, MULTIPLICATIVE, UNKNOWN };

 private:
  std::vector<double> _weights{};
  Type _type{};

 public:
  Weights() = default;
  explicit Weights(const BinTableLazy &bins, const Dataset &dset, bool rescale = false);
  Weights(const BinTableLazy &bins, const Dataset &dset, Type type, bool rescale = false);
  Weights(std::vector<double> weights, Type type) noexcept;
  Weights(std::vector<double> weights, std::string_view name);

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] double operator[](std::size_t i) const noexcept;

  [[nodiscard]] double at(std::size_t i) const;

  [[nodiscard]] const std::vector<double> &operator()() const noexcept;
  [[nodiscard]] constexpr auto type() const noexcept -> Type;

  [[nodiscard]] static auto infer_type(std::string_view name) -> Type;
  [[nodiscard]] static auto infer_type(const Dataset &dset) -> Type;
};

template <typename N>
class Balancer {
 public:
  class iterator;

 private:
  PixelSelector<N> _selector;
  std::shared_ptr<Weights> _weights;

 public:
  Balancer() = delete;
  Balancer(PixelSelector<N> selector, std::shared_ptr<Weights> weights);

  [[nodiscard]] Weights::Type type() const noexcept;

  [[nodiscard]] auto begin() const -> iterator;
  [[nodiscard]] auto end() const -> iterator;

  [[nodiscard]] auto cbegin() const -> iterator;
  [[nodiscard]] auto cend() const -> iterator;

  class iterator {
    typename PixelSelector<N>::iterator _it{};
    std::shared_ptr<Weights> _weights{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Pixel<double>;
    using pointer = value_type *;
    using reference = value_type &;
    // using iterator_category = std::random_access_iterator_tag;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    iterator(typename PixelSelector<N>::iterator it, std::shared_ptr<Weights> weights);

    [[nodiscard]] constexpr bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] constexpr bool operator<(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator<=(const iterator &other) const noexcept;

    [[nodiscard]] constexpr bool operator>(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator>=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> value_type;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
    // auto operator+=(std::size_t i) -> iterator &;
    // [[nodiscard]] auto operator+(std::size_t i) const -> iterator;

    // auto operator--() -> iterator &;
    // auto operator--(int) -> iterator;
    // auto operator-=(std::size_t i) -> iterator &;
    // [[nodiscard]] auto operator-(std::size_t i) const -> iterator;
    // [[nodiscard]] auto operator-(const iterator &other) const -> difference_type;
  };
};

}  // namespace coolerpp

#include "../../balancing_impl.hpp"
