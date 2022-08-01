#pragma once

#include <cstdint>
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

  PixelCoordinates _coords;
  const Index *_index{};
  const Dataset *_pixels_bin1_id{};
  const Dataset *_pixels_bin2_id{};
  const Dataset *_pixels_count{};
  bool _empty{true};

 public:
  PixelSelector() = delete;
  PixelSelector(const Index &index, const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                const Dataset &pixels_count, PixelCoordinates coords);

  [[nodiscard]] bool constexpr empty() const noexcept;

  [[nodiscard]] auto begin() const -> iterator;
  [[nodiscard]] auto end() const -> iterator;

  [[nodiscard]] auto cbegin() const -> iterator;
  [[nodiscard]] auto cend() const -> iterator;

  [[nodiscard]] constexpr const PixelCoordinates &coords() const noexcept;

  [[nodiscard]] static PixelCoordinates parse_query(const BinTableLazy &bins,
                                                    std::string_view query);

 private:
  [[nodiscard]] static std::uint64_t compute_end_offset(const Index &index,
                                                        const Dataset &bin2_id_dset,
                                                        std::uint64_t last_bin2_id);

 public:
  class iterator {
    friend PixelSelector<N>;
    const Index *_index{};

    std::uint32_t _first_chrom_id{};
    std::uint32_t _last_chrom_id{};

    std::uint64_t _first_bin_id{};
    std::uint64_t _last_bin_id{};

    Dataset::iterator<std::uint64_t> _bin1_id_it{};
    Dataset::iterator<std::uint64_t> _bin2_id_it{};
    Dataset::iterator<std::uint64_t> _bin2_id_last{};
    Dataset::iterator<N> _count_it{};

    explicit iterator(const Index &index, const Dataset &pixels_bin1_id,
                      const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                      const PixelCoordinates &coords) noexcept;

    [[nodiscard]] static auto make_end_iterator(const Index &index, const Dataset &pixels_bin2_id,
                                                const PixelCoordinates &coords) -> iterator;

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
  };
};

}  // namespace coolerpp

#include "../../pixel_selector_impl.hpp"
