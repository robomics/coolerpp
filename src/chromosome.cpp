// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/chromosome.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "coolerpp/common.hpp"

namespace coolerpp {

Chromosome::Chromosome(std::string name_, std::uint32_t size_) noexcept
    : name(std::move(name_)), size(size_) {}

bool Chromosome::operator<(const Chromosome& other) const noexcept {
  return this->name < other.name;
}

bool Chromosome::operator>(const Chromosome& other) const noexcept {
  return this->name > other.name;
}

bool Chromosome::operator<=(const Chromosome& other) const noexcept {
  return this->name <= other.name;
}

bool Chromosome::operator>=(const Chromosome& other) const noexcept {
  return this->name >= other.name;
}

bool Chromosome::operator==(const Chromosome& other) const noexcept {
  // Comparing by ID only is not enough, as different sets of chromosomes may use the same ID for
  // different chromosomes
  return this->name == other.name && this->size == other.size;
}

bool Chromosome::operator!=(const Chromosome& other) const noexcept { return !(*this == other); }

bool operator<(const Chromosome& a, std::string_view b_name) noexcept { return a.name < b_name; }
bool operator>(const Chromosome& a, std::string_view b_name) noexcept { return a.name > b_name; }
bool operator<=(const Chromosome& a, std::string_view b_name) noexcept { return a.name <= b_name; }
bool operator>=(const Chromosome& a, std::string_view b_name) noexcept { return a.name >= b_name; }
bool operator==(const Chromosome& a, std::string_view b_name) noexcept { return a.name == b_name; }
bool operator!=(const Chromosome& a, std::string_view b_name) noexcept { return a.name != b_name; }

bool operator<(std::string_view a_name, const Chromosome& b) noexcept { return a_name < b.name; }
bool operator>(std::string_view a_name, const Chromosome& b) noexcept { return a_name > b.name; }
bool operator<=(std::string_view a_name, const Chromosome& b) noexcept { return a_name <= b.name; }
bool operator>=(std::string_view a_name, const Chromosome& b) noexcept { return a_name >= b.name; }
bool operator==(std::string_view a_name, const Chromosome& b) noexcept { return a_name == b.name; }
bool operator!=(std::string_view a_name, const Chromosome& b) noexcept { return a_name != b.name; }

ChromosomeSet::ChromosomeSet(std::initializer_list<Chromosome> chromosomes)
    : ChromosomeSet(chromosomes.begin(), chromosomes.end()) {}

auto ChromosomeSet::begin() const -> const_iterator { return this->cbegin(); }
auto ChromosomeSet::end() const -> const_iterator { return this->cend(); }
auto ChromosomeSet::cbegin() const -> const_iterator { return this->_set.cbegin(); }
auto ChromosomeSet::cend() const -> const_iterator { return this->_set.cend(); }

auto ChromosomeSet::rbegin() const -> const_reverse_iterator { return this->rcbegin(); }
auto ChromosomeSet::rend() const -> const_reverse_iterator { return this->rcend(); }
auto ChromosomeSet::rcbegin() const -> const_reverse_iterator { return this->_set.rcbegin(); }
auto ChromosomeSet::rcend() const -> const_reverse_iterator { return this->_set.rcend(); }

bool ChromosomeSet::empty() const noexcept { return this->size() == 0; }
std::size_t ChromosomeSet::size() const noexcept { return this->_set.size(); }

auto ChromosomeSet::find(std::uint32_t id) const -> const_iterator {
  if (static_cast<std::size_t>(id) > this->size()) {
    return this->end();
  }
  return this->_set.begin() + static_cast<std::ptrdiff_t>(id);
}

auto ChromosomeSet::find(std::string_view chrom_name) const -> const_iterator {
  return this->_set.find(chrom_name);
}

auto ChromosomeSet::find(const Chromosome& chrom) const -> const_iterator {
  return this->_set.find(chrom);
}

const Chromosome& ChromosomeSet::at(std::uint32_t id) const {
  this->validate_chrom_id(id);
  return *this->find(id);
}

const Chromosome& ChromosomeSet::at(std::string_view chrom_name) const {
  if (const auto match = this->find(chrom_name); match != this->end()) {
    return *match;
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
}

bool ChromosomeSet::contains(std::uint32_t id) const { return this->find(id) != this->end(); }
bool ChromosomeSet::contains(const Chromosome& chrom) const {
  return this->find(chrom) != this->end();
}
bool ChromosomeSet::contains(std::string_view chrom_name) const {
  return this->find(chrom_name) != this->end();
}

std::uint32_t ChromosomeSet::get_id(const Chromosome& chrom) const {
  if (const auto match = this->find(chrom); match != this->end()) {
    return static_cast<std::uint32_t>(std::distance(this->begin(), match));
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome {} not found"), chrom));
}

std::uint32_t ChromosomeSet::get_id(std::string_view chrom_name) const {
  if (const auto match = this->find(chrom_name); match != this->end()) {
    return static_cast<std::uint32_t>(std::distance(this->begin(), match));
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
}

bool ChromosomeSet::operator==(const ChromosomeSet& other) const {
  return this->_set == other._set;
}

bool ChromosomeSet::operator!=(const ChromosomeSet& other) const { return !(*this == other); }

void ChromosomeSet::validate_chrom_id(std::uint32_t chrom_id) const {
  if (static_cast<std::size_t>(chrom_id) >= this->size()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome with id {} not found"), chrom_id));
  }
}

namespace internal {

bool ChromEq::operator()(const Chromosome& a, const Chromosome& b) const noexcept {
  return a.name == b.name;
}
bool ChromEq::operator()(const Chromosome& a, std::string_view b_name) const noexcept {
  return a.name == b_name;
}
bool ChromEq::operator()(std::string_view a_name, const Chromosome& b) const noexcept {
  return a_name == b.name;
}

std::size_t ChromHasher::operator()(const Chromosome& c) const {
  return std::hash<Chromosome>{}(c);
}
std::size_t ChromHasher::operator()(std::string_view name) const {
  return std::hash<std::string_view>{}(name);
}
}  // namespace internal

}  // namespace coolerpp

std::size_t std::hash<coolerpp::Chromosome>::operator()(const coolerpp::Chromosome& k) const {
  return std::hash<std::string>{}(k.name);
}
