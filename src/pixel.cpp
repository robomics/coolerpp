// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/pixel.hpp"

#include <cassert>
#include <cstdint>
#include <string_view>

#include "coolerpp/bin_table.hpp"
#include "coolerpp/chromosome.hpp"

namespace coolerpp {

PixelCoordinates::PixelCoordinates(std::shared_ptr<const BinTableLazy> bins,
                                   const Chromosome &chrom1, const Chromosome &chrom2,
                                   std::uint32_t bin1_start_, std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins->chromosomes().get_id(chrom1), bins->chromosomes().get_id(chrom2),
                       bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(std::shared_ptr<const BinTableLazy> bins,
                                   std::string_view chrom1_name, std::string_view chrom2_name,
                                   std::uint32_t bin1_start_, std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins->chromosomes().get_id(chrom1_name),
                       bins->chromosomes().get_id(chrom2_name), bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(std::shared_ptr<const BinTableLazy> bins,
                                   std::uint32_t chrom1_id_, std::uint32_t chrom2_id_,
                                   std::uint32_t bin1_start_, std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins->coord_to_bin_id(chrom1_id_, bin1_start_),
                       bins->coord_to_bin_id(chrom2_id_, bin2_start_)) {}

PixelCoordinates::PixelCoordinates(std::shared_ptr<const BinTableLazy> bins,
                                   const Chromosome &chrom, std::uint32_t bin1_start_,
                                   std::uint32_t bin2_start_)
    : PixelCoordinates(std::move(bins), chrom, chrom, bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(std::shared_ptr<const BinTableLazy> bins, std::uint32_t chrom_id,
                                   std::uint32_t bin1_start_, std::uint32_t bin2_start_)
    : PixelCoordinates(std::move(bins), chrom_id, chrom_id, bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(std::shared_ptr<const BinTableLazy> bins,
                                   std::string_view chrom_name, std::uint32_t bin1_start_,
                                   std::uint32_t bin2_start_)
    : PixelCoordinates(std::move(bins), chrom_name, chrom_name, bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(std::shared_ptr<const BinTableLazy> bins, std::uint64_t bin1_id_,
                                   std::uint64_t bin2_id_)
    : _bins(std::move(bins)), _bin1_id(bin1_id_), _bin2_id(bin2_id_) {
  assert(_bin1_id <= _bins->size());
  assert(_bin2_id <= _bins->size());
}

const Chromosome &PixelCoordinates::chrom1() const { return this->bin1().chrom; }

const Chromosome &PixelCoordinates::chrom2() const { return this->bin2().chrom; }

std::uint32_t PixelCoordinates::chrom1_id() const {
  return this->_bins->chromosomes().get_id(this->chrom1());
}

std::uint32_t PixelCoordinates::chrom2_id() const {
  return this->_bins->chromosomes().get_id(this->chrom2());
}

Bin PixelCoordinates::bin1() const {
  assert(this->_bins);
  assert(!!*this);

  return this->_bins->bin_id_to_coords(_bin1_id);
}

Bin PixelCoordinates::bin2() const {
  assert(this->_bins);
  assert(!!*this);

  return this->_bins->bin_id_to_coords(_bin2_id);
}

std::uint64_t PixelCoordinates::bin1_id() const noexcept { return this->_bin1_id; }
std::uint64_t PixelCoordinates::bin2_id() const noexcept { return this->_bin2_id; }

std::uint32_t PixelCoordinates::bin_size() const noexcept {
  assert(this->_bins);
  return this->_bins->bin_size();
}

}  // namespace coolerpp
