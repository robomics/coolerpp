// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <variant>

#include "coolerpp/common.hpp"
#include "group.hpp"

namespace coolerpp {

struct Attribute {
  // Variants are listed in order from the most common to the least common for perf. reasons
  // clang-format off
  using AttributeVar = std::variant<
      std::monostate,
      std::string,
      std::uint64_t,
      std::int64_t,
      double,
      std::uint32_t, std::uint16_t, std::uint8_t,
      std::int32_t, std::int16_t, std::int8_t,
      float>;
  // clang-format on
  Attribute() = delete;

  // ParentObj e.g. DataSet, Group
  template <class ParentObj>
  [[nodiscard]] static bool exists(ParentObj& h5obj, std::string_view key);

  template <class T, class ParentObj>
  static void write(ParentObj& h5obj, std::string_view key, const T& value,
                    bool overwrite_if_exists = false);
  template <class ParentObj>
  [[nodiscard]] static auto read(const ParentObj& h5obj, std::string_view key,
                                 bool missing_ok = false) -> AttributeVar;
  template <class T, class ParentObj>
  [[nodiscard]] static T read(const ParentObj& h5obj, std::string_view key);

  template <class T, class ParentObj>
  [[nodiscard]] static std::vector<T> read_vector(const ParentObj& h5obj, std::string_view key);
  template <class T, class ParentObj>
  static void read_vector(const ParentObj& h5obj, std::string_view key, std::vector<T>& buff);

 private:
  template <std::size_t i = 1>  // i = 1 skips T=monostate
  [[nodiscard]] static auto read_variant(const HighFive::Attribute& attr) -> AttributeVar;
  template <class T1, class Tout, class Tin = remove_cvref_t<T1>>
  [[nodiscard]] static Tout numeric_converter(T1& buff);
};
}  // namespace coolerpp

#include "../../attribute_impl.hpp"
