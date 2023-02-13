// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <memory>
#include <type_traits>

#include "coolerpp/common.hpp"
#include "coolerpp/dataset.hpp"
#include "coolerpp/pixel.hpp"

namespace coolerpp {

class BinTableLazy;
class Index;

template <class N>
class PixelSelector {
  static_assert(std::is_arithmetic_v<N>);

 public:
  class iterator;

 private:
  static constexpr std::size_t chunk_size = 4096ULL << 10U;

  std::shared_ptr<PixelCoordinates> _coord1{};
  std::shared_ptr<PixelCoordinates> _coord2{};
  const Index *_index{};
  const Dataset *_pixels_bin1_id{};
  const Dataset *_pixels_bin2_id{};
  const Dataset *_pixels_count{};

  PixelSelector(const Index &index, const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                const Dataset &pixels_count, std::shared_ptr<PixelCoordinates> coords) noexcept;
  PixelSelector(const Index &index, const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                const Dataset &pixels_count, std::shared_ptr<PixelCoordinates> coord1,
                std::shared_ptr<PixelCoordinates> coord2) noexcept;

 public:
  PixelSelector() = delete;
  PixelSelector(const Index &index, const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                const Dataset &pixels_count) noexcept;
  PixelSelector(const Index &index, const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                const Dataset &pixels_count, PixelCoordinates coords) noexcept;

  PixelSelector(const Index &index, const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                const Dataset &pixels_count, PixelCoordinates coord1,
                PixelCoordinates coord2) noexcept;

  [[nodiscard]] auto begin() const -> iterator;
  [[nodiscard]] auto end() const -> iterator;

  [[nodiscard]] auto cbegin() const -> iterator;
  [[nodiscard]] auto cend() const -> iterator;

  [[nodiscard]] constexpr const PixelCoordinates &coord1() const noexcept;
  [[nodiscard]] constexpr const PixelCoordinates &coord2() const noexcept;

  [[nodiscard]] static PixelCoordinates parse_query(const BinTableLazy &bins,
                                                    std::string_view query);

  class iterator {
    friend PixelSelector<N>;
    const Index *_index{};

    std::shared_ptr<PixelCoordinates> _coord1{};
    std::shared_ptr<PixelCoordinates> _coord2{};

    Dataset::iterator<std::uint64_t> _bin1_id_it{};
    Dataset::iterator<std::uint64_t> _bin2_id_it{};
    Dataset::iterator<std::uint64_t> _bin2_id_last{};
    Dataset::iterator<N> _count_it{};

    explicit iterator(const Index &index, const Dataset &pixels_bin1_id,
                      const Dataset &pixels_bin2_id, const Dataset &pixels_count);

    explicit iterator(const Index &index, const Dataset &pixels_bin1_id,
                      const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                      std::shared_ptr<PixelCoordinates> coord1,
                      std::shared_ptr<PixelCoordinates> coord2);

    static auto at_end(const Index &index, const Dataset &pixels_bin1_id,
                       const Dataset &pixels_bin2_id, const Dataset &pixels_count) -> iterator;

    static auto at_end(const Index &index, const Dataset &pixels_bin1_id,
                       const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                       std::shared_ptr<PixelCoordinates> coord1,
                       std::shared_ptr<PixelCoordinates> coord2) -> iterator;

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Pixel<N>;
    using pointer = value_type *;
    using reference = value_type &;
    // using iterator_category = std::random_access_iterator_tag;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;

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
   private:
    void jump_to_row(std::uint64_t bin_id);
    void jump_to_col(std::uint64_t bin_id);
    void jump(std::uint64_t bin1_id, std::uint64_t bin2_id);
  };
};

}  // namespace coolerpp

#include "../../pixel_selector_impl.hpp"
