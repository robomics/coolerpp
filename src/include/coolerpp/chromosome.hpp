// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <tsl/ordered_set.h>

#include <cstdint>
#include <initializer_list>
#include <string>
#include <vector>

namespace coolerpp {
struct Chromosome;
}

namespace std {
template <>
struct hash<coolerpp::Chromosome> {
  size_t operator()(const coolerpp::Chromosome& k) const;
};
}  // namespace std

namespace fmt {
template <>
struct formatter<coolerpp::Chromosome> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  template <typename FormatContext>
  auto format(const coolerpp::Chromosome& c, FormatContext& ctx) const -> decltype(ctx.out());
};
}  // namespace fmt

namespace coolerpp {

namespace internal {
struct ChromEq {
  using is_transparent = void;
  bool operator()(const Chromosome& a, const Chromosome& b) const noexcept;
  bool operator()(const Chromosome& a, std::string_view b_name) const noexcept;
  bool operator()(std::string_view a_name, const Chromosome& b) const noexcept;
};

struct ChromHasher {
  std::size_t operator()(const Chromosome& c) const;
  std::size_t operator()(std::string_view name) const;
};
}  // namespace internal

struct Chromosome {
  std::string name{};
  std::uint32_t size{};

  Chromosome() = default;
  Chromosome(std::string name_, std::uint32_t size_) noexcept;

  [[nodiscard]] bool operator<(const Chromosome& other) const noexcept;
  [[nodiscard]] bool operator>(const Chromosome& other) const noexcept;
  [[nodiscard]] bool operator<=(const Chromosome& other) const noexcept;
  [[nodiscard]] bool operator>=(const Chromosome& other) const noexcept;
  [[nodiscard]] bool operator==(const Chromosome& other) const noexcept;
  [[nodiscard]] bool operator!=(const Chromosome& other) const noexcept;

  friend bool operator<(const Chromosome& a, std::string_view b_name) noexcept;
  friend bool operator>(const Chromosome& a, std::string_view b_name) noexcept;
  friend bool operator<=(const Chromosome& a, std::string_view b_name) noexcept;
  friend bool operator>=(const Chromosome& a, std::string_view b_name) noexcept;
  friend bool operator==(const Chromosome& a, std::string_view b_name) noexcept;
  friend bool operator!=(const Chromosome& a, std::string_view b_name) noexcept;

  friend bool operator<(std::string_view a_name, const Chromosome& b) noexcept;
  friend bool operator>(std::string_view a_name, const Chromosome& b) noexcept;
  friend bool operator<=(std::string_view a_name, const Chromosome& b) noexcept;
  friend bool operator>=(std::string_view a_name, const Chromosome& b) noexcept;
  friend bool operator==(std::string_view a_name, const Chromosome& b) noexcept;
  friend bool operator!=(std::string_view a_name, const Chromosome& b) noexcept;
};

class ChromosomeSet {
  using SetT = tsl::ordered_set<Chromosome, internal::ChromHasher, internal::ChromEq>;
  SetT _set{};

 public:
  using key_type = typename SetT::key_type;
  using value_type = typename SetT::value_type;
  using size_type = typename SetT::size_type;
  using difference_type = typename SetT::difference_type;
  using hasher = typename SetT::hasher;
  using key_equal = typename SetT::key_equal;
  using allocator_type = typename SetT::allocator_type;
  using reference = typename SetT::const_reference;
  using const_reference = typename SetT::const_reference;
  using pointer = typename SetT::const_pointer;
  using const_pointer = typename SetT::const_pointer;
  using iterator = typename SetT::const_iterator;
  using const_iterator = typename SetT::const_iterator;
  using reverse_iterator = typename SetT::reverse_iterator;
  using const_reverse_iterator = typename SetT::const_reverse_iterator;

  ChromosomeSet() = default;

  template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
  ChromosomeSet(ChromosomeNameIt first_chrom_name, ChromosomeNameIt last_chrom_name,
                ChromosomeSizeIt first_chrom_size);

  // Note: chromosome IDs are not preserved
  template <typename ChromosomeIt>
  ChromosomeSet(ChromosomeIt first_chrom, ChromosomeIt last_chrom);
  ChromosomeSet(std::initializer_list<Chromosome> chromosomes);

  [[nodiscard]] auto begin() const -> const_iterator;
  [[nodiscard]] auto end() const -> const_iterator;
  [[nodiscard]] auto cbegin() const -> const_iterator;
  [[nodiscard]] auto cend() const -> const_iterator;

  [[nodiscard]] auto rbegin() const -> const_reverse_iterator;
  [[nodiscard]] auto rend() const -> const_reverse_iterator;
  [[nodiscard]] auto rcbegin() const -> const_reverse_iterator;
  [[nodiscard]] auto rcend() const -> const_reverse_iterator;

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;

  [[nodiscard]] auto find(std::uint32_t id) const -> const_iterator;
  [[nodiscard]] auto find(std::string_view chrom_name) const -> const_iterator;
  [[nodiscard]] auto find(const Chromosome& chrom) const -> const_iterator;

  [[nodiscard]] const Chromosome& at(std::uint32_t id) const;
  [[nodiscard]] const Chromosome& at(std::string_view chrom_name) const;

  [[nodiscard]] const Chromosome& operator[](std::uint32_t id) const noexcept;
  [[nodiscard]] const Chromosome& operator[](std::string_view chrom_name) const noexcept;

  [[nodiscard]] bool contains(std::uint32_t id) const;
  [[nodiscard]] bool contains(const Chromosome& chrom) const;
  [[nodiscard]] bool contains(std::string_view chrom_name) const;

  [[nodiscard]] std::uint32_t get_id(const Chromosome& chrom) const;
  [[nodiscard]] std::uint32_t get_id(std::string_view chrom_name) const;

  [[nodiscard]] bool operator==(const ChromosomeSet& other) const;
  [[nodiscard]] bool operator!=(const ChromosomeSet& other) const;

  [[nodiscard]] const Chromosome& find_longest_chromosome() const;
  [[nodiscard]] const Chromosome& find_chromosome_with_longest_name() const;

 private:
  void validate_chrom_id(std::uint32_t chrom_id) const;

  template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
  [[nodiscard]] static auto construct_set(ChromosomeNameIt first_chrom_name,
                                          ChromosomeNameIt last_chrom_name,
                                          ChromosomeSizeIt first_chrom_size) -> SetT;

  template <typename ChromosomeIt>
  [[nodiscard]] static auto construct_set(ChromosomeIt first_chrom, ChromosomeIt last_chrom)
      -> SetT;
};
}  // namespace coolerpp

#include "../../chromosome_impl.hpp"
