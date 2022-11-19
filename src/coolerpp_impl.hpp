// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5Utility.hpp>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

#include "coolerpp/attribute.hpp"
#include "coolerpp/bin_table.hpp"
#include "coolerpp/dataset.hpp"
#include "coolerpp/group.hpp"
#include "coolerpp/internal/numeric_utils.hpp"
#include "coolerpp/internal/prime_number_table.hpp"
#include "coolerpp/internal/type_pretty_printer.hpp"
#include "coolerpp/internal/variant_buff.hpp"
#include "coolerpp/uri.hpp"
#include "coolerpp/utils.hpp"

namespace coolerpp {

template <class InputIt>
inline void init_mcool(std::string_view file_path, InputIt first_resolution,
                       InputIt last_resolution, bool force_overwrite) {
  using I = remove_cvref_t<decltype(*first_resolution)>;
  static_assert(std::is_integral_v<I>,
                "InputIt should be an iterator over a collection of integral numbers.");
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  const auto mode = force_overwrite ? HighFive::File::Truncate : HighFive::File::Create;
  HighFive::File fp(std::string{file_path}, mode);
  Attribute::write(fp, "format", std::string{MCOOL_MAGIC});
  Attribute::write(fp, "format-version", std::int64_t(3));

  auto res_group = fp.createGroup("/resolutions");
  std::for_each(first_resolution, last_resolution, [&](auto res) {
    assert(res > 0);
    auto cooler_root = res_group.createGroup(fmt::to_string(res));
  });
}

// template <class ChromSizeInputIt, class CellIDInputIt>
// inline void init_scool(std::string_view file_path, ChromSizeInputIt first_chrom,
//                        ChromSizeInputIt last_chrom, CellIDInputIt first_cell_id,
//                        CellIDInputIt last_cell_id, std::uint32_t bin_size, bool force_overwrite)
//                        {
//   [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
//   const auto mode = force_overwrite ? IO_MODE::Truncate : IO_MODE::Create;
//   HighFive::File fp(std::string{file_path}, static_cast<unsigned>(mode));
//   fp.createAttribute("format", std::string{SCOOL_MAGIC});
//   fp.createAttribute("format-version", std::int64_t(3));
//
//   const auto bin_table = binnify(first_chrom, last_chrom, bin_size);
//
//   auto res_group = fp.createGroup("/resolutions");
//
// }

template <class PixelT>
inline File::File(std::string_view uri, ChromosomeSet chroms, [[maybe_unused]] PixelT pixel,
                  StandardAttributes attributes)
    : _mode(HighFive::File::ReadWrite),
      _fp(std::make_unique<HighFive::File>(open_file(uri, _mode, false))),
      _root_group(open_or_create_root_group(*_fp, uri)),
      _groups(create_groups(_root_group)),
      _datasets(create_datasets<PixelT>(_root_group, chroms)),
      _attrs(std::move(attributes)),
      _pixel_variant(PixelT(0)),
      _bins(std::make_unique<BinTable>(std::move(chroms), this->bin_size())),
      _index(std::make_unique<Index>(*_bins)) {
  assert(this->bin_size() != 0);
  assert(!_bins->empty());
  assert(!chromosomes().empty());
  assert(!_index->empty());
  assert(std::holds_alternative<PixelT>(this->_pixel_variant));

  this->write_sentinel_attr();
}

template <class PixelT>
inline File File::create_new_cooler(std::string_view uri, const ChromosomeSet &chroms,
                                    std::uint32_t bin_size, bool overwrite_if_exists,
                                    StandardAttributes attributes) {
  static_assert(std::is_arithmetic_v<PixelT>);
  assert(bin_size != 0);
  attributes.bin_size = bin_size;
  try {
    const auto [file_path, root_path] = parse_cooler_uri(uri);
    const auto uri_is_file_path = root_path.empty() || root_path == "/";

    // URI is like myfile.mcool::/resolutions/100, but myfile.mcool does not exist
    if (!uri_is_file_path && !std::filesystem::exists(file_path)) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("parent file \"{}\" does not exist.\n"
                     "Did you forget to create the parent file with e.g. init_mcool()?"),
          uri, file_path));
    }

    // URI points to an existing file, but overwrite_if_exists=false
    if (!overwrite_if_exists && uri_is_file_path && std::filesystem::exists(file_path)) {
      throw std::runtime_error("URI points to an existing file");
    }

    auto mode = overwrite_if_exists ? HighFive::File::Overwrite : HighFive::File::Create;

    // File exists but cooler may not
    if (std::filesystem::exists(file_path) && !uri_is_file_path) {
      mode = HighFive::File::ReadWrite;
    }

    {
      auto fp = open_file(uri, mode, false);
      auto root_group = open_or_create_root_group(fp, uri);
      if (!uri_is_file_path && utils::is_cooler(root_group())) {
        if (overwrite_if_exists) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("overwriting cooler nested inside .mcool or .scool is not yet supported.\n"
                         "Path to parent file: \"{}\"\""
                         "Path to nested cooler: \"{}\""),
              file_path, root_path));
        }
      }
      assert(!utils::is_cooler(root_group()));
    }
    // At this point the parent file is guaranteed to exist, so we can always open it in ReadWrite
    // mode
    return File(uri, chroms, PixelT(0), attributes);

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Cannot create cooler at the following URI: \"{}\". Reason: {}"),
                    uri, e.what()));
  }
}

template <class PixelT>
inline void File::create(std::string_view uri, const coolerpp::ChromosomeSet &chroms,
                         std::uint32_t bin_size, bool overwrite_if_exists,
                         coolerpp::StandardAttributes attributes) {
  *this = File::create_new_cooler<PixelT>(uri, chroms, bin_size, overwrite_if_exists, attributes);
}

namespace internal {
template <class Variant, std::size_t i = 0>
[[nodiscard]] inline Variant read_pixel_variant(const HighFive::DataSet &dset) {
  if constexpr (i < std::variant_size_v<Variant>) {
    using T = std::variant_alternative_t<i, Variant>;
    if (dset.getDataType() != HighFive::create_datatype<T>()) {
      return read_pixel_variant<Variant, i + 1>(dset);
    }
    return T{};
  }

  constexpr bool variant_has_monostate =
      std::is_same_v<std::monostate, std::variant_alternative_t<0, Variant>>;
  if constexpr (variant_has_monostate) {
    return std::monostate();
  } else {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unsupported type for dataset \"{}\""), dset.getPath()));
  }
}
}  // namespace internal

template <class PixelT>
inline auto File::create_datasets(RootGroup &root_grp, const ChromosomeSet &chroms) -> DatasetMap {
  DatasetMap datasets(MANDATORY_DATASET_NAMES.size() + 1);

  auto create_dataset = [&](const auto &path, const auto &type) {
    using T = remove_cvref_t<decltype(type)>;
    if constexpr (is_string_v<T>) {
      const auto chrom_with_longest_name = find_chromosome_with_longest_name(chroms);
      datasets.emplace(path, Dataset{root_grp, path, chrom_with_longest_name.name,
                                     HighFive::DataSpace::UNLIMITED});
    } else {
      datasets.emplace(path, Dataset{root_grp, path, type, HighFive::DataSpace::UNLIMITED});
    }
  };

  create_dataset("chroms/name", std::string{});
  create_dataset("chroms/length", std::int32_t{});

  create_dataset("bins/chrom", std::int32_t{});
  create_dataset("bins/start", std::int32_t{});
  create_dataset("bins/end", std::int32_t{});

  create_dataset("pixels/bin1_id", std::int64_t{});
  create_dataset("pixels/bin2_id", std::int64_t{});
  create_dataset("pixels/count", PixelT{});

  create_dataset("indexes/bin1_offset", std::int64_t{});
  create_dataset("indexes/chrom_offset", std::int64_t{});

  assert(datasets.size() == MANDATORY_DATASET_NAMES.size());

  return datasets;
}

template <class PixelIt>
inline void File::validate_pixels_before_append(PixelIt first_pixel, PixelIt last_pixel) const {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using T = decltype(std::declval<PixelT>().count);
  try {
    std::for_each(first_pixel, last_pixel, [&](const Pixel<T> &pixel) {
      if (pixel.count == T{0}) {
        throw std::runtime_error("found a pixel of value 0");
      }

      if (!this->chromosomes().contains(*pixel.coords.chrom1)) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("invalid chromosome id {}"), *pixel.coords.chrom1));
      }

      if (*pixel.coords.chrom1 != *pixel.coords.chrom2 &&
          !this->chromosomes().contains(*pixel.coords.chrom2)) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("invalid chromosome id {}"), *pixel.coords.chrom2));
      }

      if (const auto bin_id = pixel.coords.bin1_id(); bin_id > this->bins().size()) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("invalid bin id {}: bin maps outside of the bin table"), bin_id));
      }

      if (const auto bin_id = pixel.coords.bin2_id(); bin_id > this->bins().size()) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("invalid bin id {}: bin maps outside of the bin table"), bin_id));
      }
    });

    if (!this->dataset("pixels/bin1_id").empty()) {
      const auto last_bin1 = this->dataset("pixels/bin1_id").read_last<std::size_t>();
      const auto last_bin2 = this->dataset("pixels/bin2_id").read_last<std::size_t>();

      const auto new_bin1 = first_pixel->coords.bin1_id();
      const auto new_bin2 = first_pixel->coords.bin2_id();

      if (last_bin1 == new_bin1) {
        if (last_bin2 >= new_bin2) {
          const auto coord1 = this->bins().bin_id_to_coords(new_bin2);
          const auto coord2 = this->bins().bin_id_to_coords(last_bin2);
          throw std::runtime_error(fmt::format(
              FMT_STRING("new pixel {} is located upstream of pixel {}"), coord1, coord2));
        }
      } else if (last_bin1 >= new_bin1) {
        const auto coord1 = this->bins().bin_id_to_coords(new_bin1);
        const auto coord2 = this->bins().bin_id_to_coords(last_bin1);
        throw std::runtime_error(fmt::format(
            FMT_STRING("new pixel {} is located upstream of pixel {}"), coord1, coord2));
      }
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("pixel validation failed: {}"), e.what()));
  }
}

template <class T>
inline bool File::has_pixel_of_type() const noexcept {
  return std::holds_alternative<T>(this->_pixel_variant);
}

template <class PixelT, class>
inline StandardAttributes StandardAttributes::init(std::uint32_t bin_size_) {
  StandardAttributes attrs{};
  attrs.bin_size = bin_size_;
  if constexpr (std::is_floating_point_v<PixelT>) {
    attrs.sum = 0.0;
  } else if constexpr (std::is_signed_v<PixelT>) {
    attrs.sum = std::int64_t(0);
  } else {
    assert(std::is_unsigned_v<PixelT>);
    attrs.sum = std::uint64_t(0);
  }
  return attrs;
}

template <class ChromIt, class UnaryOperation, class>
inline void File::write_chromosomes(Dataset &name_dset, Dataset &size_dset, ChromIt first_chrom,
                                    ChromIt last_chrom, UnaryOperation op) {
  const auto num_chroms = std::distance(first_chrom, last_chrom);
  if (num_chroms == 0) {
    return;
  }

  try {
    name_dset.write(first_chrom, last_chrom, 0, true,
                    [&](const auto &chrom) { return op(chrom).name; });
  } catch (const HighFive::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to write {} chromosome name(s) to \"{}\": {}"), num_chroms,
                    name_dset.uri(), e.what()));
  }
  try {
    size_dset.write(first_chrom, last_chrom, 0, true,
                    [&](const auto &chrom) { return op(chrom).size; });
  } catch (const HighFive::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to write {} chromosome size(s) to \"{}\": {}"), num_chroms,
                    size_dset.uri(), e.what()));
  }

  assert(name_dset.size() == static_cast<std::size_t>(num_chroms));
  assert(size_dset.size() == static_cast<std::size_t>(num_chroms));
}

template <class PixelIt>
inline void File::update_indexes(PixelIt first_pixel, PixelIt last_pixel) {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using T = decltype(std::declval<PixelT>().count);

  if (first_pixel == last_pixel) {
    return;
  }

  const auto last_bin_written = this->get_last_bin_written();

  auto nnz = *this->_attrs.nnz;
  PixelCoordinates first_pixel_in_row{this->bins(), last_bin_written.chrom.name,
                                      last_bin_written.bin_start, last_bin_written.bin_start};

  std::for_each(first_pixel, last_pixel, [&](const Pixel<T> &p) {
    if (first_pixel_in_row.bin1_start != p.coords.bin1_start) {
      first_pixel_in_row = p.coords;
      this->index().set_offset_by_bin_id(first_pixel_in_row.bin1_id(), nnz);
    }
    nnz++;
  });
}

template <class PixelIt, class>
inline void File::append_pixels(PixelIt first_pixel, PixelIt last_pixel, bool validate,
                                [[maybe_unused]] std::size_t chunk_size) {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using T = decltype(std::declval<PixelT>().count);

  this->update_indexes(first_pixel, last_pixel);

  if (validate) {
    this->validate_pixels_before_append(first_pixel, last_pixel);
  }

  this->dataset("pixels/bin1_id").append(first_pixel, last_pixel, [&](const Pixel<T> &pixel) {
    return pixel.coords.bin1_id();
  });

  this->dataset("pixels/bin2_id").append(first_pixel, last_pixel, [&](const Pixel<T> &pixel) {
    return pixel.coords.bin2_id();
  });

  T sum = 0;
  this->dataset("pixels/count").append(first_pixel, last_pixel, [&](const Pixel<T> &pixel) {
    sum += pixel.count;
    return pixel.count;
  });

  this->_attrs.nnz = this->dataset("pixels/bin1_id").size();

  if constexpr (std::is_floating_point_v<T>) {
    std::get<double>(this->_attrs.sum) += sum;
  } else if constexpr (std::is_signed_v<T>) {
    std::get<std::int64_t>(this->_attrs.sum) += sum;
  } else {
    assert(std::is_unsigned_v<T>);
    std::get<std::uint64_t>(this->_attrs.sum) += sum;
  }
}

template <class N>
inline typename PixelSelector<N>::iterator File::begin() const {
  return this
      ->fetch<N>(this->chromosomes().begin()->name, 0, this->chromosomes().rbegin()->name,
                 this->chromosomes().rbegin()->size - 1)
      .begin();
}

template <class N>
inline typename PixelSelector<N>::iterator File::cbegin() const {
  return this->begin<N>();
}

template <class N>
inline typename PixelSelector<N>::iterator File::end() const {
  return this
      ->fetch<N>(this->chromosomes().begin()->name, 0, this->chromosomes().rbegin()->name,
                 this->chromosomes().rbegin()->size - 1)
      .end();
}

template <class N>
inline typename PixelSelector<N>::iterator File::cend() const {
  return this->end<N>();
}

template <class N>
inline PixelSelector<N> File::fetch(std::string_view query) const {
  // clang-format off
  return PixelSelector<N>(this->index(),
                          this->dataset("pixels/bin1_id"),
                          this->dataset("pixels/bin2_id"),
                          this->dataset("pixels/count"),
                          PixelSelector<N>::parse_query(this->bins(), query));
  // clang-format on
}

template <class N>
inline PixelSelector<N> File::fetch(std::string_view chrom1_name, std::uint32_t pos1,
                                    std::string_view chrom2_name, std::uint32_t pos2) const {
  // clang-format off
  return PixelSelector<N>(this->index(),
                          this->dataset("pixels/bin1_id"),
                          this->dataset("pixels/bin2_id"),
                          this->dataset("pixels/count"),
                          PixelCoordinates{this->bins(), chrom1_name, chrom2_name, pos1, pos2});
  // clang-format on
}

}  // namespace coolerpp
