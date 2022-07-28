// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <fmt/ranges.h>

namespace coolerpp::utils {

constexpr ValidationStatusCooler::operator bool() const noexcept { return this->is_cooler; }

constexpr ValidationStatusMultiresCooler::operator bool() const noexcept {
  return this->is_multires_file;
}

constexpr ValidationStatusScool::operator bool() const noexcept { return this->is_scool_file; }

}  // namespace coolerpp::utils

constexpr auto fmt::formatter<coolerpp::utils::ValidationStatusCooler>::parse(
    format_parse_context &ctx) -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.end();
}

template <class FormatContext>
auto fmt::formatter<coolerpp::utils::ValidationStatusCooler>::format(
    const coolerpp::utils::ValidationStatusCooler &s, FormatContext &ctx) const
    -> decltype(ctx.out()) {
  auto bool_to_str = [](bool x) { return x ? "true" : "false"; };

  // clang-format off
  return fmt::format_to(
      ctx.out(),
      FMT_STRING("uri=\"{}\"\n"
                 "is_hdf5={}\n"
                 "missing_or_invalid_format_attr={}\n"
                 "missing_or_invalid_bin_type_attr={}\n"
                 "missing_groups=[{}]\n"
                 "is_valid_cooler={}"),
      s.uri,
      bool_to_str(s.is_hdf5),
      bool_to_str(s.missing_or_invalid_format_attr),
      bool_to_str(s.missing_or_invalid_bin_type_attr),
      fmt::join(s.missing_groups, ", "),
      bool_to_str(s.is_cooler));
  // clang-format on
}

constexpr auto fmt::formatter<coolerpp::utils::ValidationStatusMultiresCooler>::parse(
    format_parse_context &ctx) -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.end();
}

template <class FormatContext>
auto fmt::formatter<coolerpp::utils::ValidationStatusMultiresCooler>::format(
    const coolerpp::utils::ValidationStatusMultiresCooler &s, FormatContext &ctx) const
    -> decltype(ctx.out()) {
  auto bool_to_str = [](bool x) { return x ? "true" : "false"; };

  // clang-format off
  return fmt::format_to(
      ctx.out(),
      FMT_STRING("uri=\"{}\"\n"
                 "is_hdf5={}\n"
                 "missing_or_invalid_format_attr={}\n"
                 "missing_or_invalid_bin_type_attr={}\n"
                 "missing_groups=[{}]\n"
                 "is_valid_multires_file={}\n"
                 "invalid_resolutions{}{}"),
      s.uri,
      bool_to_str(s.is_hdf5),
      bool_to_str(s.missing_or_invalid_format_attr),
      bool_to_str(s.missing_or_invalid_bin_type_attr),
      fmt::join(s.missing_groups, ", "),
      bool_to_str(s.is_multires_file),
      s.invalid_resolutions.empty() ? "=[]" : ":\n - ",
      fmt::join(s.invalid_resolutions, "\n - "));
  // clang-format on
}

constexpr auto fmt::formatter<coolerpp::utils::ValidationStatusScool>::parse(
    format_parse_context &ctx) -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.end();
}

template <class FormatContext>
auto fmt::formatter<coolerpp::utils::ValidationStatusScool>::format(
    const coolerpp::utils::ValidationStatusScool &s, FormatContext &ctx) const
    -> decltype(ctx.out()) {
  auto bool_to_str = [](bool x) { return x ? "true" : "false"; };

  // clang-format off
  return fmt::format_to(
      ctx.out(),
      FMT_STRING("uri=\"{}\"\n"
                 "is_hdf5={}\n"
                 "missing_or_invalid_format_attr={}\n"
                 "missing_or_invalid_bin_type_attr={}\n"
                 "missing_groups=[{}]\n"
                 "is_valid_scool_file={}\n"
                 "invalid_cells{}{}"),
      s.uri,
      bool_to_str(s.is_hdf5),
      bool_to_str(s.missing_or_invalid_format_attr),
      bool_to_str(s.missing_or_invalid_bin_type_attr),
      fmt::join(s.missing_groups, ", "),
      bool_to_str(s.is_scool_file),
      s.invalid_cells.empty() ? "=[]" : ":\n - ",
      fmt::join(s.invalid_cells, "\n - "));
  // clang-format on
}
