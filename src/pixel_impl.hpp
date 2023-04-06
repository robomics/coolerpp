// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cassert>
#include <cstdint>
#include <string_view>

#include "coolerpp/bin_table.hpp"
#include "coolerpp/chromosome.hpp"

namespace coolerpp {

inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTableLazy> &bins,
                                          const Chromosome &chrom1, const Chromosome &chrom2,
                                          std::uint32_t bin1_start_, std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins->chromosomes().get_id(chrom1), bins->chromosomes().get_id(chrom2),
                       bin1_start_, bin2_start_) {}

inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTableLazy> &bins,
                                          std::string_view chrom1_name,
                                          std::string_view chrom2_name, std::uint32_t bin1_start_,
                                          std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins->chromosomes().get_id(chrom1_name),
                       bins->chromosomes().get_id(chrom2_name), bin1_start_, bin2_start_) {}

inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTableLazy> &bins,
                                          std::uint32_t chrom1_id_, std::uint32_t chrom2_id_,
                                          std::uint32_t bin1_start_, std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins->coord_to_bin_id(chrom1_id_, bin1_start_),
                       bins->coord_to_bin_id(chrom2_id_, bin2_start_)) {}

inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTableLazy> &bins,
                                          const Chromosome &chrom, std::uint32_t bin1_start_,
                                          std::uint32_t bin2_start_)
    : PixelCoordinates(bins, chrom, chrom, bin1_start_, bin2_start_) {}

inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTableLazy> &bins,
                                          std::uint32_t chrom_id, std::uint32_t bin1_start_,
                                          std::uint32_t bin2_start_)
    : PixelCoordinates(bins, chrom_id, chrom_id, bin1_start_, bin2_start_) {}

inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTableLazy> &bins,
                                          std::string_view chrom_name, std::uint32_t bin1_start_,
                                          std::uint32_t bin2_start_)
    : PixelCoordinates(bins, chrom_name, chrom_name, bin1_start_, bin2_start_) {}

inline PixelCoordinates::PixelCoordinates(std::shared_ptr<const BinTableLazy> bins,
                                          std::uint64_t bin1_id_, std::uint64_t bin2_id_)
    : _bins(std::move(bins)), _bin1_id(bin1_id_), _bin2_id(bin2_id_) {
  assert(_bin1_id <= _bins->size());
  assert(_bin2_id <= _bins->size());
}

inline PixelCoordinates::operator bool() const noexcept { return !!this->_bins; }

inline const Chromosome &PixelCoordinates::chrom1() const { return this->bin1().chrom; }

inline const Chromosome &PixelCoordinates::chrom2() const { return this->bin2().chrom; }

inline std::uint32_t PixelCoordinates::chrom1_id() const {
  return this->_bins->chromosomes().get_id(this->chrom1());
}

inline std::uint32_t PixelCoordinates::chrom2_id() const {
  return this->_bins->chromosomes().get_id(this->chrom2());
}

inline Bin PixelCoordinates::bin1() const {
  assert(this->_bins);
  assert(!!*this);

  return this->_bins->bin_id_to_coords(_bin1_id);
}

inline Bin PixelCoordinates::bin2() const {
  assert(this->_bins);
  assert(!!*this);

  return this->_bins->bin_id_to_coords(_bin2_id);
}

inline std::uint64_t PixelCoordinates::bin1_id() const noexcept { return this->_bin1_id; }
inline std::uint64_t PixelCoordinates::bin2_id() const noexcept { return this->_bin2_id; }

inline std::uint32_t PixelCoordinates::bin_size() const noexcept {
  assert(this->_bins);
  return this->_bins->bin_size();
}

constexpr bool PixelCoordinates::operator==(const PixelCoordinates &other) const noexcept {
  return this->_bin1_id == other._bin1_id && this->_bin2_id == other._bin2_id;
}

constexpr bool PixelCoordinates::operator!=(const PixelCoordinates &other) const noexcept {
  return !(*this == other);
}

constexpr bool PixelCoordinates::operator<(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id < other._bin2_id;
  }
  return this->_bin1_id < other._bin1_id;
}

constexpr bool PixelCoordinates::operator<=(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id <= other._bin2_id;
  }
  return this->_bin1_id <= other._bin1_id;
}

constexpr bool PixelCoordinates::operator>(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id > other._bin2_id;
  }
  return this->_bin1_id > other._bin1_id;
}

constexpr bool PixelCoordinates::operator>=(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id >= other._bin2_id;
  }
  return this->_bin1_id >= other._bin1_id;
}

template <typename N>
inline Pixel<N>::operator bool() const noexcept {
  return !!this->coords;
}
template <typename N>
constexpr bool Pixel<N>::operator==(const Pixel<N> &other) const noexcept {
  return this->coords == other.coords && this->count == other.count;
}
template <typename N>
constexpr bool Pixel<N>::operator!=(const Pixel<N> &other) const noexcept {
  return !(*this == other);
}
template <typename N>
constexpr bool Pixel<N>::operator<(const Pixel<N> &other) const noexcept {
  return this->coords < other.coords;
}
template <typename N>
constexpr bool Pixel<N>::operator<=(const Pixel<N> &other) const noexcept {
  return this->coords <= other.coords;
}
template <typename N>
constexpr bool Pixel<N>::operator>(const Pixel<N> &other) const noexcept {
  return this->coords > other.coords;
}
template <typename N>
constexpr bool Pixel<N>::operator>=(const Pixel<N> &other) const noexcept {
  return this->coords >= other.coords;
}

}  // namespace coolerpp

constexpr auto fmt::formatter<coolerpp::PixelCoordinates>::parse(format_parse_context &ctx)
    -> decltype(ctx.begin()) {
  const auto *it = ctx.begin();
  const auto *end = ctx.end();
  const auto fmt_string =
      std::string_view{&(*ctx.begin()), static_cast<std::size_t>(ctx.end() - ctx.begin())};

  if (it != end) {
    if (fmt_string.find("bedpe") != std::string_view::npos) {
      this->presentation = Presentation::bedpe;
      it += std::string_view{"bedpe"}.size();  // NOLINT
    } else if (fmt_string.find("raw") != std::string_view::npos) {
      this->presentation = Presentation::raw;
      it += std::string_view{"raw"}.size();  // NOLINT
    }
  }

  if (it != end && *it != '}') {
    throw fmt::format_error("invalid format");
  }

  return it;
}

template <typename FormatContext>
inline auto fmt::formatter<coolerpp::PixelCoordinates>::format(const coolerpp::PixelCoordinates &c,
                                                               FormatContext &ctx) const
    -> decltype(ctx.out()) {
  if (this->presentation == Presentation::raw) {
    // clang-format off
    return fmt::format_to(ctx.out(),
                          FMT_STRING("{}\t{}"),
                          c.bin1_id(), c.bin2_id());
    // clang-format on
  }

  assert(this->presentation == Presentation::bedpe);
  const auto bin1 = c.bin1();
  const auto bin2 = c.bin2();
  // clang-format off
  return fmt::format_to(ctx.out(),
                        FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}"),
                        bin1.chrom.name,
                        bin1.start,
                        (std::min)(bin1.start + c.bin_size(), bin1.chrom.size),
                        bin2.chrom.name,
                        bin2.start,
                        (std::min)(bin2.start + c.bin_size(), bin2.chrom.size));
  // clang-format on
}

template <typename N>
constexpr auto fmt::formatter<coolerpp::Pixel<N>>::parse(format_parse_context &ctx)
    -> decltype(ctx.begin()) {
  this->coord_formatter.parse(ctx);
  return ctx.end();
}

template <typename N>
constexpr auto fmt::formatter<coolerpp::Pixel<N>>::presentation() const noexcept -> Presentation {
  return this->coord_formatter.presentation;
}

template <typename N>
template <typename FormatContext>
inline auto fmt::formatter<coolerpp::Pixel<N>>::format(const coolerpp::Pixel<N> &p,
                                                       FormatContext &ctx) const
    -> decltype(ctx.out()) {
  if (this->presentation() == Presentation::raw) {
    return fmt::format_to(ctx.out(), FMT_STRING("{:raw}\t{}"), p.coords, p.count);
  }

  assert(this->presentation() == Presentation::bedpe);
  return fmt::format_to(ctx.out(), FMT_STRING("{:bedpe}\t{}"), p.coords, p.count);
}
