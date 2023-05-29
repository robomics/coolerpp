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

template <typename N, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
class PixelSelector {
  static_assert(std::is_arithmetic_v<N>);

 public:
  class iterator;

 private:
  std::shared_ptr<PixelCoordinates> _coord1{};
  std::shared_ptr<PixelCoordinates> _coord2{};
  std::shared_ptr<const Index> _index{};
  const Dataset *_pixels_bin1_id{};
  const Dataset *_pixels_bin2_id{};
  const Dataset *_pixels_count{};

  PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                const std::shared_ptr<PixelCoordinates> &coords) noexcept;
  PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                std::shared_ptr<PixelCoordinates> coord1,
                std::shared_ptr<PixelCoordinates> coord2) noexcept;

 public:
  PixelSelector() = delete;
  PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                const Dataset &pixels_bin2_id, const Dataset &pixels_count) noexcept;
  PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                PixelCoordinates coords) noexcept;

  PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                const Dataset &pixels_bin2_id, const Dataset &pixels_count, PixelCoordinates coord1,
                PixelCoordinates coord2) noexcept;

  template <std::size_t CHUNK_SIZE_OTHER>
  [[nodiscard]] bool operator==(const PixelSelector<N, CHUNK_SIZE_OTHER> &other) const noexcept;
  template <std::size_t CHUNK_SIZE_OTHER>
  [[nodiscard]] bool operator!=(const PixelSelector<N, CHUNK_SIZE_OTHER> &other) const noexcept;

  [[nodiscard]] auto begin() const -> iterator;
  [[nodiscard]] auto end() const -> iterator;

  [[nodiscard]] auto cbegin() const -> iterator;
  [[nodiscard]] auto cend() const -> iterator;

  [[nodiscard]] constexpr const PixelCoordinates &coord1() const noexcept;
  [[nodiscard]] constexpr const PixelCoordinates &coord2() const noexcept;

  [[nodiscard]] static PixelCoordinates parse_query(std::shared_ptr<const BinTableLazy> bins,
                                                    std::string_view query);

  class iterator {
    friend PixelSelector<N, CHUNK_SIZE>;
    std::shared_ptr<const Index> _index{};

    std::shared_ptr<PixelCoordinates> _coord1{};
    std::shared_ptr<PixelCoordinates> _coord2{};

    Dataset::iterator<std::uint64_t, CHUNK_SIZE> _bin1_id_it{};
    Dataset::iterator<std::uint64_t, CHUNK_SIZE> _bin2_id_it{};
    Dataset::iterator<std::uint64_t, CHUNK_SIZE> _bin2_id_last{};
    Dataset::iterator<N, CHUNK_SIZE> _count_it{};

    explicit iterator(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                      const Dataset &pixels_bin2_id, const Dataset &pixels_count);

    explicit iterator(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                      const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                      std::shared_ptr<PixelCoordinates> coord1,
                      std::shared_ptr<PixelCoordinates> coord2);

    static auto at_end(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                       const Dataset &pixels_bin2_id, const Dataset &pixels_count) -> iterator;

    static auto at_end(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
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
    void jump_to_next_overlap();
    [[nodiscard]] bool discard() const;
    [[nodiscard]] std::size_t h5_offset() const noexcept;
    void jump_at_end();
    void refresh();
  };
};

}  // namespace coolerpp

#include "../../pixel_selector_impl.hpp"
