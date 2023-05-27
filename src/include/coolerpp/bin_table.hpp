// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "coolerpp/chromosome.hpp"
#include "coolerpp/common.hpp"

namespace coolerpp {

struct Bin {
  const Chromosome &chrom;
  std::uint32_t start;
  std::uint32_t end;

  Bin() = delete;
  constexpr Bin(const Chromosome &chrom_, std::uint32_t start_, std::uint32_t end_) noexcept;
  [[nodiscard]] bool operator==(const Bin &other) const noexcept;
  [[nodiscard]] bool operator!=(const Bin &other) const noexcept;

  [[nodiscard]] bool operator<(const Bin &other) const noexcept;
  [[nodiscard]] bool operator<=(const Bin &other) const noexcept;

  [[nodiscard]] bool operator>(const Bin &other) const noexcept;
  [[nodiscard]] bool operator>=(const Bin &other) const noexcept;
};

struct BinTable {
  std::vector<const Chromosome *> chroms{};
  std::vector<std::uint32_t> bin_starts{};
  std::vector<std::uint32_t> bin_ends{};
};

class BinTableLazy {
  ChromosomeSet _chroms{};
  std::vector<std::uint64_t> _num_bins_prefix_sum{};
  std::uint32_t _bin_size{std::numeric_limits<std::uint32_t>::max()};

 public:
  class iterator;
  using const_iterator = const iterator;
  friend iterator;

  BinTableLazy() = default;
  BinTableLazy(ChromosomeSet chroms, std::uint32_t bin_size);
  template <typename ChromIt>
  BinTableLazy(ChromIt first_chrom, ChromIt last_chrom, std::uint32_t bin_size);
  template <typename ChromNameIt, typename ChromSizeIt>
  BinTableLazy(ChromNameIt first_chrom_name, ChromNameIt last_chrom_name,
               ChromSizeIt first_chrom_size, std::uint32_t bin_size);

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t num_chromosomes() const;
  [[nodiscard]] constexpr std::uint32_t bin_size() const noexcept;
  [[nodiscard]] constexpr const ChromosomeSet &chromosomes() const noexcept;

  [[nodiscard]] constexpr const std::vector<std::uint64_t> &num_bin_prefix_sum() const noexcept;

  [[nodiscard]] constexpr auto begin() -> iterator;
  [[nodiscard]] constexpr auto end() -> iterator;
  [[nodiscard]] constexpr auto begin() const -> const_iterator;
  [[nodiscard]] constexpr auto end() const -> const_iterator;
  [[nodiscard]] constexpr auto cbegin() const -> const_iterator;
  [[nodiscard]] constexpr auto cend() const -> const_iterator;

  [[nodiscard]] BinTableLazy at(const Chromosome &chrom) const;
  [[nodiscard]] BinTableLazy at(std::string_view chrom_name) const;
  [[nodiscard]] BinTableLazy at(std::uint32_t chrom_id) const;

  // Map bin_id to chromosome + (relative) position
  [[nodiscard]] Bin bin_id_to_coords(std::uint64_t bin_id) const;

  // Map genomic coords to bin_id
  [[nodiscard]] std::uint64_t coord_to_bin_id(const Bin &bin) const;
  [[nodiscard]] std::uint64_t coord_to_bin_id(const Chromosome &chrom, std::uint32_t pos) const;
  [[nodiscard]] std::uint64_t coord_to_bin_id(std::string_view chrom_name, std::uint32_t pos) const;
  [[nodiscard]] std::uint64_t coord_to_bin_id(std::uint32_t chrom_id, std::uint32_t pos) const;

  [[nodiscard]] BinTable concretize() const;

  [[nodiscard]] bool operator==(const BinTableLazy &other) const;
  [[nodiscard]] bool operator!=(const BinTableLazy &other) const;

 private:
  [[nodiscard]] static std::vector<std::uint64_t> compute_num_bins_prefix_sum(
      const ChromosomeSet &chroms, std::uint32_t bin_size);

 public:
  class iterator {
    friend BinTableLazy;
    const BinTableLazy *_bin_table{};
    std::size_t _idx{0};
    std::uint32_t _chrom_id{0};

    static constexpr auto npos = std::numeric_limits<std::size_t>::max();

    constexpr explicit iterator(const BinTableLazy &bin_table) noexcept;

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Bin;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::bidirectional_iterator_tag;

    constexpr iterator() noexcept = default;

    constexpr bool operator==(const iterator &other) const noexcept;
    constexpr bool operator!=(const iterator &other) const noexcept;
    auto operator*() const -> value_type;
    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
    auto operator--() -> iterator &;
    auto operator--(int) -> iterator;

   private:
    [[nodiscard]] static constexpr auto make_end_iterator(const BinTableLazy &table) noexcept
        -> iterator;
    [[nodiscard]] const Chromosome &chromosome(std::uint32_t chrom_id) const;
    [[nodiscard]] const Chromosome &chromosome() const;
    [[nodiscard]] constexpr std::uint32_t bin_size() const noexcept;
    [[nodiscard]] std::uint64_t compute_num_bins() const noexcept;
    [[nodiscard]] std::size_t num_chromosomes() const noexcept;
  };
};

}  // namespace coolerpp

namespace fmt {
template <>
struct formatter<coolerpp::Bin> {
  constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin());
  template <typename FormatContext>
  auto format(const coolerpp::Bin &b, FormatContext &ctx) const -> decltype(ctx.out());
};
}  // namespace fmt

#include "../../bin_table_impl.hpp"
