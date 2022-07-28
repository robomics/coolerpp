// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <tsl/hopscotch_map.h>

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <string>

#include "coolerpp/internal/suppress_warnings.hpp"

namespace coolerpp {

DISABLE_WARNING_PUSH
DISABLE_WARNING_DEPRECATED_DECLARATIONS
struct RootGroup {
  HighFive::Group group{};

  [[nodiscard]] constexpr HighFive::Group &operator()() noexcept { return this->group; };
  [[nodiscard]] constexpr const HighFive::Group &operator()() const noexcept {
    return this->group;
  };

  [[nodiscard]] inline std::string file_name() const { return this->group.getFile().getName(); }
  [[nodiscard]] inline std::string hdf5_path() const { return this->group.getPath(); }
  [[nodiscard]] inline std::string uri() const {
    return fmt::format(FMT_STRING("{}::{}"), this->file_name(), this->hdf5_path());
  }
};

struct Group {
  RootGroup root_group{};
  HighFive::Group group{};

  [[nodiscard]] constexpr HighFive::Group &operator()() noexcept { return this->group; };
  [[nodiscard]] constexpr const HighFive::Group &operator()() const noexcept {
    return this->group;
  };
};
DISABLE_WARNING_POP

using GroupMap = tsl::hopscotch_map<std::string, Group>;

}  // namespace coolerpp
