// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <highfive/H5Exception.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "coolerpp/attribute.hpp"
#include "coolerpp/bin_table.hpp"
#include "coolerpp/chromosome.hpp"
#include "coolerpp/dataset.hpp"
#include "coolerpp/group.hpp"
#include "coolerpp/internal/type_pretty_printer.hpp"
#include "coolerpp/uri.hpp"

namespace coolerpp {

inline bool File::check_sentinel_attr(const HighFive::Group &grp) {
  const auto generated_by_v = Attribute::read(grp, "generated-by", true);
  if (const auto *generated_by = std::get_if<std::string>(&generated_by_v);
      !generated_by || generated_by->find("coolerpp") == std::string::npos) {
    return false;
  }

  using T = remove_cvref_t<decltype(internal::SENTINEL_ATTR_VALUE)>;
  const auto sentinel_v = Attribute::read(grp, internal::SENTINEL_ATTR_NAME, true);
  const auto *sentinel = std::get_if<T>(&sentinel_v);

  return static_cast<bool>(sentinel) && *sentinel == internal::SENTINEL_ATTR_VALUE;
}

namespace internal {
[[nodiscard]] inline std::vector<std::uint64_t> import_chrom_offsets(const Dataset &dset,
                                                                     std::size_t expected_size) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto offsets = dset.read_all<std::vector<std::uint64_t>>();
  try {
    if (offsets.size() != expected_size) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("expected {} offsets, found {}"), expected_size, offsets.size()));
    }
    if (offsets.front() != 0) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("first offset should be 0, found {}"), offsets.front()));
    }
    if (!std::is_sorted(offsets.begin(), offsets.end())) {
      throw std::runtime_error("offsets are not in ascending order");
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to import offsets from {}: {}"), dset.uri(), e.what()));
  }

  return offsets;
}
}  // namespace internal

template <class N>
inline PixelSelector<N> File::fetch(std::string_view query) const {
  return this->fetch<N>(PixelSelector<N>::parse_query(this->_bins, query));
}

template <class N>
inline PixelSelector<N> File::fetch(std::string_view chrom, std::uint32_t start,
                                    std::uint32_t end) const {
  return this->fetch<N>(PixelCoordinates{this->_bins, chrom, start, end - std::min(1U, end)});
}

template <class N>
inline PixelSelector<N> File::fetch(PixelCoordinates query) const {
  // clang-format off
  return PixelSelector<N>(this->_index,
                          this->dataset("pixels/bin1_id"),
                          this->dataset("pixels/bin2_id"),
                          this->dataset("pixels/count"),
                          std::move(query));
  // clang-format on
}

template <class N>
inline PixelSelector<N> File::fetch(std::string_view range1, std::string_view range2) const {
  if (range1 == range2) {
    return this->fetch<N>(range1);
  }

  return this->fetch<N>(PixelSelector<N>::parse_query(this->_bins, range1),
                        PixelSelector<N>::parse_query(this->_bins, range2));
}

template <class N>
inline PixelSelector<N> File::fetch(std::string_view chrom1, std::uint32_t start1,
                                    std::uint32_t end1, std::string_view chrom2,
                                    std::uint32_t start2, std::uint32_t end2) const {
  // clang-format off
  return PixelSelector<N>(this->_index,
                          this->dataset("pixels/bin1_id"),
                          this->dataset("pixels/bin2_id"),
                          this->dataset("pixels/count"),
                          PixelCoordinates{this->_bins, chrom1, start1, end1},
                          PixelCoordinates{this->_bins, chrom2, start2, end2});
  // clang-format on
}

template <class N>
inline PixelSelector<N> File::fetch(PixelCoordinates coord1, PixelCoordinates coord2) const {
  // clang-format off
  return PixelSelector<N>(this->_index,
                          this->dataset("pixels/bin1_id"),
                          this->dataset("pixels/bin2_id"),
                          this->dataset("pixels/count"),
                          std::move(coord1), std::move(coord2));
  // clang-format on
}

inline std::shared_ptr<Weights> File::read_weights(std::string_view name) {
  if (name.empty()) {
    throw std::runtime_error("weight dataset name is empty");
  }

  return this->read_weights(name, Weights::infer_type(name));
}

inline std::shared_ptr<Weights> File::read_weights(std::string_view name, Weights::Type type) {
  if (name.empty()) {
    throw std::runtime_error("weight dataset name is empty");
  }

  const auto dset_path =
      fmt::format(FMT_STRING("{}/{}"), this->_groups.at("bins").group.getPath(), name);
  if (const auto it = this->_weights.find(dset_path); it != this->_weights.end()) {
    return it->second;
  }

  if (!this->_root_group().exist(dset_path)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to read \"{}\" weights: dataset \"{}\" does not exist"),
                    name, dset_path));
  }

  const auto node = this->_weights.emplace(
      name, std::make_shared<Weights>(*this->_bins, Dataset{this->_root_group, dset_path}, type));
  return node.first->second;
}

inline bool File::purge_weights(std::string_view name) {
  if (this->_weights.empty()) {
    return false;
  }
  if (name == "") {
    this->_weights.clear();
    return true;
  }
  return this->_weights.erase(std::string{name});
}

inline auto File::open_root_group(const HighFive::File &f, std::string_view uri) -> RootGroup {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  return {f.getGroup(parse_cooler_uri(uri).group_path)};
}

inline auto File::open_groups(const RootGroup &root_grp) -> GroupMap {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  GroupMap groups(MANDATORY_GROUP_NAMES.size() + 1);
  groups.emplace(root_grp.hdf5_path(), Group{root_grp, root_grp()});

  std::transform(MANDATORY_GROUP_NAMES.begin(), MANDATORY_GROUP_NAMES.end(),
                 std::inserter(groups, groups.begin()), [&root_grp](const auto group_name) {
                   const auto name = std::string{group_name};
                   auto group_obj = root_grp().getGroup(std::string{group_name});

                   return std::make_pair(name, Group{root_grp, group_obj});
                 });
  return groups;
}

inline auto File::open_datasets(const RootGroup &root_grp, std::string_view weight_dataset)
    -> DatasetMap {
  DatasetMap datasets(MANDATORY_DATASET_NAMES.size() + 1);

  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto open_dataset = [&](const auto dataset_uri) {
    return std::make_pair(std::string{dataset_uri}, Dataset{root_grp, dataset_uri});
  };

  std::transform(MANDATORY_DATASET_NAMES.begin(), MANDATORY_DATASET_NAMES.end(),
                 std::inserter(datasets, datasets.begin()), open_dataset);

  const auto path = fmt::format(FMT_STRING("bins/{}"), weight_dataset);
  if (root_grp().exist(path)) {
    datasets.emplace(open_dataset(path));
  }

  return datasets;
}

inline auto File::read_standard_attributes(const RootGroup &root_grp, bool initialize_missing)
    -> StandardAttributes {
  auto attrs = initialize_missing ? StandardAttributes::init(0) : StandardAttributes::init_empty();
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT

  auto read_or_throw = [&](const auto &key, auto &buff) {
    using T = remove_cvref_t<decltype(buff)>;
    try {
      buff = Attribute::read<T>(root_grp(), key);
    } catch (const std::exception &e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Failed to read attribute \"{}\" from path \"{}\". Reason: {}"),
                      key, root_grp().getPath(), e.what()));
    }
  };

  auto read_optional = [&](const auto &key, auto &buff, bool missing_ok) {
    if (!Attribute::exists(root_grp(), key) && missing_ok) {
      return false;
    }

    try {
      using T = remove_cvref_t<decltype(*buff)>;
      buff = Attribute::read<T>(root_grp(), key);
      return true;
    } catch (const std::exception &e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Failed to read attribute \"{}\" from path \"{}\". Reason: {}"),
                      key, root_grp().getPath(), e.what()));
    }
  };

  auto read_sum_optional = [&](bool missing_ok, std::string_view key, auto &buff) {
    if (!Attribute::exists(root_grp(), key) && missing_ok) {
      return false;
    }

    try {
      auto sumv = Attribute::read(root_grp(), key);
      std::visit(
          [&](auto sum) {
            using T = remove_cvref_t<decltype(sum)>;
            if constexpr (std::is_integral_v<T>) {
              buff = conditional_static_cast<std::int64_t>(sum);
              return;
            }
            if constexpr (std::is_floating_point_v<T>) {
              buff = conditional_static_cast<double>(sum);
              return;
            }
            throw std::runtime_error(
                fmt::format(FMT_STRING("Attribute \"{}{}\" as an unexpected type. Expected a "
                                       "numeric type, found {}"),
                            root_grp().getPath(), key, internal::type_name<T>()));
          },
          sumv);
      return true;
    } catch (const std::exception &e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Failed to read attribute \"{}\" from path \"{}\". Reason: {}"),
                      key, root_grp().getPath(), e.what()));
    }
  };

  // Read mandatory attributes
  // We read format-version first because some attributes are mandatory only for cooler v3
  read_or_throw("format-version", attrs.format_version);
  read_or_throw("bin-size", attrs.bin_size);
  read_or_throw("format", attrs.format);

  // Read mandatory attributes for Cooler v3
  auto missing_ok = attrs.format_version < 3;
  read_optional("bin-type", attrs.bin_type, missing_ok);
  read_optional("storage-mode", attrs.storage_mode, missing_ok);

  // Try to read reserved attributes
  missing_ok = true;
  read_optional("creation-date", attrs.creation_date, missing_ok);
  read_optional("format-url", attrs.format_url, missing_ok);
  read_optional("generated-by", attrs.generated_by, missing_ok);

  if (!read_optional("genome-assembly", attrs.assembly, missing_ok)) {
    read_optional("assembly", attrs.assembly, missing_ok);
  }

  read_optional("metadata", attrs.metadata, missing_ok);

  // Try to read other common attributes
  read_optional("nbins", attrs.nbins, missing_ok);
  read_optional("nchroms", attrs.nchroms, missing_ok);
  read_optional("nnz", attrs.nnz, missing_ok);

  read_sum_optional(missing_ok, "sum", attrs.sum);
  read_sum_optional(missing_ok, "cis", attrs.cis);

  return attrs;
}

inline auto File::import_chroms(const Dataset &chrom_names, const Dataset &chrom_sizes,
                                bool missing_ok) -> ChromosomeSet {
  try {
    [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
    std::vector<std::string> names;
    std::vector<std::uint32_t> sizes;
    chrom_names.read_all(names);
    chrom_sizes.read_all(sizes);

    if (names.size() != sizes.size()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Cooler file \"{}\" appears to be corrupted: {} and "
                                 "{} shape mismatch: found {} name(s) and {} length(s)"),
                      chrom_names.file_name(), chrom_names.hdf5_path(), chrom_sizes.hdf5_path(),
                      names.size(), sizes.size()));
    }

    return ChromosomeSet{names.begin(), names.end(), sizes.begin()};

  } catch ([[maybe_unused]] const HighFive::Exception &e) {
    if (missing_ok) {
      return {};
    }
    throw;
  }
}

inline Index File::import_indexes(const Dataset &chrom_offset_dset, const Dataset &bin_offset_dset,
                                  const ChromosomeSet &chroms,
                                  std::shared_ptr<const BinTableLazy> bin_table,
                                  std::uint64_t expected_nnz, bool missing_ok) {
  assert(bin_table);
  try {
    if (bin_offset_dset.empty()) {
      assert(chrom_offset_dset.empty());
      if (missing_ok) {
        return Index{bin_table};
      }
      throw std::runtime_error("index datasets are empty");
    }

    if (bin_offset_dset.size() != bin_table->size() + 1) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("failed to import offsets from {}: expected {} offsets, found {}"),
                      bin_offset_dset.hdf5_path(), bin_table->size() + 1, bin_offset_dset.size()));
    }

    const auto chrom_offsets = internal::import_chrom_offsets(chrom_offset_dset, chroms.size() + 1);

    Index idx{bin_table, expected_nnz};

    std::size_t bin_id = 0;
    std::for_each(bin_offset_dset.begin<std::uint64_t>(), bin_offset_dset.end<std::uint64_t>(),
                  [&](std::uint64_t offset) {
                    if (bin_id < bin_table->size()) {
                      idx.set_offset_by_bin_id(bin_id++, offset);
                    } else {
                      // Last bin
                      assert(bin_id == bin_table->size());
                    }
                  });

    try {
      idx.validate();
    } catch (const std::exception &e) {
      throw std::runtime_error(fmt::format(FMT_STRING("index validation failed: {}"), e.what()));
    }

    return idx;

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to import indexes for cooler at URI: \"{}\": {}"),
                    bin_offset_dset.get_parent().uri(), e.what()));
  }
}

inline bool File::check_sentinel_attr() { return File::check_sentinel_attr(this->_root_group()); }

inline auto File::get_last_bin_written() const -> Bin {
  const auto &dset = this->dataset("pixels/bin1_id");
  if (dset.empty()) {
    return this->bins().bin_id_to_coords(0);
  }
  const auto bin1_id = dset.read_last<std::uint64_t>();
  return this->bins().bin_id_to_coords(bin1_id);
}

}  // namespace coolerpp
