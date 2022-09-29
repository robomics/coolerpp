// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/utils.hpp"

#include <fmt/chrono.h>
#include <fmt/format.h>

#include <cstdint>
#include <highfive/H5File.hpp>
#include <highfive/H5Utility.hpp>
#include <string>
#include <string_view>

#include "coolerpp/attribute.hpp"
#include "coolerpp/common.hpp"
#include "coolerpp/uri.hpp"

namespace coolerpp::utils {

ValidationStatusCooler is_cooler(std::string_view uri) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};
  const auto [file_path, root_path] = parse_cooler_uri(uri);

  HighFive::File fp(file_path, HighFive::File::ReadOnly);
  return is_cooler(fp, root_path);
}

ValidationStatusMultiresCooler is_multires_file(std::string_view uri, bool validate_resolutions,
                                                std::int64_t min_version) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};
  const auto file_path = parse_cooler_uri(uri).file_path;

  HighFive::File fp(file_path, HighFive::File::ReadOnly);
  return is_multires_file(fp, validate_resolutions, min_version);
}

ValidationStatusScool is_scool_file(std::string_view uri, bool validate_cells) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};
  const auto file_path = parse_cooler_uri(uri).file_path;

  HighFive::File fp(file_path, HighFive::File::ReadOnly);
  return is_scool_file(fp, validate_cells);
}

ValidationStatusCooler is_cooler(const HighFive::File &fp, std::string_view root_path) {
  return is_cooler(fp.getGroup(std::string{root_path}));
}

ValidationStatusCooler is_cooler(const HighFive::Group &root_group) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  ValidationStatusCooler status{};
  status.uri = root_group.getFile().getName();
  const auto root_path = root_group.getPath();
  if (!root_path.empty() && root_path != "/") {
    // NOLINTNEXTLINE(readability-implicit-bool-conversion)
    status.uri += fmt::format(FMT_STRING("::/{}"), root_path.substr(root_path.front() == '/'));
  }

  status.is_hdf5 = root_group.getFile().isValid();
  // Check file is in HDF5 format
  if (!status.is_hdf5) {
    return status;
  }

  // Check file has the appropriate format attribute
  if (Attribute::exists(root_group, "format")) {
    const auto format = Attribute::read<std::string>(root_group, "format");
    status.missing_or_invalid_format_attr = format != COOL_MAGIC;
  }

  if (Attribute::exists(root_group, "format-version")) {
    const auto version = Attribute::read<std::uint8_t>(root_group, "format-version");
    status.file_was_properly_closed = version != coolerpp::internal::SENTINEL_ATTR_VALUE;
    status.missing_or_invalid_format_attr |= version == 0 || version > 3;
  }

  // Check file has a bin-type that we support (currently only "fixed" is supported)
  if (Attribute::exists(root_group, "bin-type")) {
    const auto bin_type = Attribute::read<std::string>(root_group, "bin-type");
    status.missing_or_invalid_bin_type_attr = bin_type != "fixed";
  }

  // Check file has the mandatory groups
  for (const auto &group : MANDATORY_GROUP_NAMES) {
    try {
      root_group.getGroup(std::string{group});
    } catch (const HighFive::Exception &e) {
      status.missing_groups.emplace_back(group);
    }
  }

  // clang-format off
  status.is_cooler = status.is_hdf5 &&
                     status.file_was_properly_closed &&
                     !status.missing_or_invalid_format_attr &&
                     !status.missing_or_invalid_bin_type_attr &&
                     status.missing_groups.empty();
  // clang-format on

  return status;
}

ValidationStatusMultiresCooler is_multires_file(const HighFive::File &fp, bool validate_resolutions,
                                                std::int64_t min_version) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  ValidationStatusMultiresCooler status{};
  status.uri = fp.getName();

  status.file_was_properly_closed = true;
  status.is_hdf5 = fp.isValid();
  // Check file is in HDF5 format
  if (!status.is_hdf5) {
    return status;
  }

  // Check file has the appropriate format attribute
  if (Attribute::exists(fp, "format")) {
    const auto format = Attribute::read<std::string>(fp, "format");
    status.missing_or_invalid_format_attr = format != MCOOL_MAGIC;
  }

  // Check file has the appropriate format attribute
  if (Attribute::exists(fp, "format-version")) {
    const auto version = Attribute::read<std::uint8_t>(fp, "format-version");
    status.missing_or_invalid_format_attr |= version == 0 || version > 3;
  }

  // Check file has a bin-type that we support (currently only "fixed" is supported)
  // NOTE: .mcool files are not required to advertise the bin type they are using at the root level
  status.missing_or_invalid_bin_type_attr = false;
  if (Attribute::exists(fp, "bin-type")) {
    const auto bin_type = Attribute::read<std::string>(fp, "bin-type");
    status.missing_or_invalid_bin_type_attr = bin_type != "fixed";
  }

  // Try to read resolutions from the Cooler's root
  const auto resolutions = [&]() -> std::vector<std::string> {
    try {
      auto root_grp = fp.getGroup("/resolutions");
      auto resolutions_ = root_grp.listObjectNames();

      // TODO: I am not sure I fully understand the docs. here:
      //       https://cooler.readthedocs.io/en/latest/schema.html#multi-resolution
      if (min_version >= 2) {
        return resolutions_;
      }

      // The old multi-resolution layout used resolutions strictly in increments of powers of 2.
      // In this layout (MCOOL version 2), the data collections are named by zoom level,
      // starting with XYZ.1000.mcool::0 being the coarsest resolution up until the finest or
      // “base” resolution (e.g., XYZ.1000.mcool::14 for 14 levels of coarsening).
      auto it = std::find(resolutions_.begin(), resolutions_.end(), "0");
      if (it != resolutions_.end()) {
        return {};
      }
      return resolutions_;
    } catch (const HighFive::GroupException &e) {
      return {};
    }
  }();

  if (resolutions.empty()) {
    status.missing_groups.emplace_back("resolutions");
  }

  if (validate_resolutions) {
    for (const auto &resolution : resolutions) {
      const auto suffix = fmt::format(FMT_STRING("resolutions/{}"), resolution);

      if (auto status_ = is_cooler(fp, suffix); !status_) {
        status.file_was_properly_closed &= status_.file_was_properly_closed;
        status.invalid_resolutions.emplace_back(std::move(status_));
      }
    }
  }

  // clang-format off
  status.is_multires_file = status.is_hdf5 &&
                            status.file_was_properly_closed &&
                            !status.missing_or_invalid_format_attr &&
                            !status.missing_or_invalid_bin_type_attr &&
                            status.missing_groups.empty() &&
                            status.invalid_resolutions.empty();
  // clang-format on

  return status;
}

ValidationStatusScool is_scool_file(const HighFive::File &fp, bool validate_cells) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  ValidationStatusScool status{};
  status.uri = fp.getName();

  status.file_was_properly_closed = true;
  status.is_hdf5 = fp.isValid();
  // Check file is in HDF5 format
  if (!status.is_hdf5) {
    return status;
  }

  // Check file has the appropriate format attribute
  if (Attribute::exists(fp, "format")) {
    const auto format = Attribute::read<std::string>(fp, "format");
    status.missing_or_invalid_format_attr = format != SCOOL_MAGIC;
  }

  if (Attribute::exists(fp, "format")) {
    const auto version = Attribute::read<std::uint8_t>(fp, "format-version");
    status.missing_or_invalid_format_attr |= version == 0 || version > 3;
  }

  // Check file has a bin-type that we support (currently only "fixed" is supported)
  // NOTE: .scool files are not required to advertise the bin type they are using at the root level
  status.missing_or_invalid_bin_type_attr = false;
  if (Attribute::exists(fp, "bin-type")) {
    const auto bin_type = Attribute::read<std::string>(fp, "bin-type");
    status.missing_or_invalid_bin_type_attr = bin_type != "fixed";
  }

  constexpr std::array<std::string_view, 3> scool_root_groups{"chroms", "bins", "cells"};

  const auto root_groups = fp.listObjectNames();

  // Check file has the mandatory groups
  for (const auto &group : scool_root_groups) {
    try {
      fp.getGroup(std::string{group});
    } catch (const HighFive::Exception &e) {
      status.missing_groups.emplace_back(group);
    }
  }

  const auto cells = [&]() -> std::vector<std::string> {
    try {
      auto root_grp = fp.getGroup("/cells");
      return root_grp.listObjectNames();
    } catch (const HighFive::GroupException &e) {
      return {};
    }
  }();

  // Try to check if the number of groups under /cells is consistent with available metadata
  status.unexpected_number_of_cells = false;
  if (Attribute::exists(fp, "ncells")) {
    const auto expected_num_cells = Attribute::read<std::uint64_t>(fp, "ncells");
    status.unexpected_number_of_cells =
        expected_num_cells != conditional_static_cast<std::uint64_t>(cells.size());
  }

  if (validate_cells) {
    for (const auto &cell : cells) {
      const auto suffix = fmt::format(FMT_STRING("cells/{}"), cell);
      if (auto status_ = is_cooler(fp, suffix); !status_) {
        status.file_was_properly_closed &= status_.file_was_properly_closed;
        status.invalid_cells.emplace_back(std::move(status_));
      }
    }
  }

  // clang-format off
  status.is_scool_file = status.is_hdf5 &&
                         status.file_was_properly_closed &&
                         !status.missing_or_invalid_format_attr &&
                         !status.missing_or_invalid_bin_type_attr &&
                         status.missing_groups.empty() &&
                         !status.unexpected_number_of_cells &&
                         status.invalid_cells.empty();
  // clang-format on

  return status;
}

}  // namespace coolerpp::utils
