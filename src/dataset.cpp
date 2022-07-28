// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/dataset.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5Exception.hpp>
#include <highfive/H5File.hpp>
#include <string>
#include <vector>

#include "coolerpp/internal/prime_number_table.hpp"

namespace coolerpp {

Dataset::Dataset(RootGroup root_group, HighFive::DataSet dset)
    : _root_group(std::move(root_group)), _dataset(std::move(dset)) {}

Dataset::Dataset(RootGroup root_group, std::string_view path_to_dataset,
                 const HighFive::DataSetAccessProps &aprops)
    : Dataset(root_group, root_group().getDataSet(std::string{path_to_dataset}, aprops)) {}

Dataset::Dataset(RootGroup root_group, std::string_view path_to_dataset,
                 std::string_view longest_str, std::size_t max_dim,
                 const HighFive::DataSetAccessProps &aprops,
                 const HighFive::DataSetCreateProps &cprops)
    : Dataset(root_group, create_fixed_str_dataset(root_group, path_to_dataset, longest_str.size(),
                                                   max_dim, aprops, cprops)) {}

std::string Dataset::file_name() const { return this->_root_group().getFile().getName(); }

std::string Dataset::hdf5_path() const { return this->_dataset.getPath(); }

std::string Dataset::uri() const {
  return fmt::format(FMT_STRING("{}::{}"), this->file_name(), this->hdf5_path());
}

std::size_t Dataset::size() const { return this->_dataset.getElementCount(); }

bool Dataset::empty() const { return this->size() == 0; }

HighFive::DataSet Dataset::get() { return this->_dataset; }
const HighFive::DataSet &Dataset::get() const { return this->_dataset; }

RootGroup Dataset::get_parent() const { return this->_root_group; }

void Dataset::resize(std::size_t new_size) {
  if (new_size > this->_dataset.getElementCount()) {
    this->_dataset.resize({new_size});
  }
}

HighFive::Selection Dataset::select(std::size_t i) {
#if defined(__GNUC__) && !defined(__clang__)
  return this->_dataset.select(std::vector<std::size_t>{i});
#else
  return this->_dataset.select({i});
#endif
}

HighFive::Selection Dataset::select(std::size_t i) const {
#if defined(__GNUC__) && !defined(__clang__)
  return this->_dataset.select(std::vector<std::size_t>{i});
#else
  return this->_dataset.select({i});
#endif
}

HighFive::Selection Dataset::select(std::size_t i1, std::size_t i2) {
#if defined(__GNUC__) && !defined(__clang__)
  return this->_dataset.select(std::vector<std::size_t>{i1}, std::vector<std::size_t>{i2});
#else
  return this->_dataset.select({i1}, {i2});
#endif
}

HighFive::Selection Dataset::select(std::size_t i1, std::size_t i2) const {
#if defined(__GNUC__) && !defined(__clang__)
  return this->_dataset.select(std::vector<std::size_t>{i1}, std::vector<std::size_t>{i2});
#else
  return this->_dataset.select({i1}, {i2});
#endif
}

std::size_t Dataset::read(std::vector<std::string> &buff, std::size_t num,
                          std::size_t offset) const {
  if (offset + num > this->size()) {
    this->throw_out_of_range_excp(offset, num);
  }

  buff.resize(num);
  const auto str_length = this->_dataset.getDataType().getSize();
  std::string strbuff(this->size() * str_length, '\0');
  this->_dataset.read(strbuff.data(), this->_dataset.getDataType());

  for (std::size_t i = 0; i < buff.size(); ++i) {
    const auto i0 = (offset + i) * str_length;
    const auto i1 = std::min(strbuff.find('\0', i0), i0 + str_length);

    buff[i] = strbuff.substr(i0, i1 - i0);
  }

  return offset + buff.size();
}

std::size_t Dataset::write(const std::vector<std::string> &buff, std::size_t offset,
                           bool allow_dataset_resize) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  if (offset + buff.size() > this->size()) {
    if (allow_dataset_resize) {
      this->resize(offset + buff.size());
    } else {
      this->throw_out_of_range_excp(offset, buff.size());
    }
  }
  // TODO construct string, then write
  const auto str_length = this->get_h5type().getSize();
  std::string strbuff(str_length * buff.size(), '\0');
  for (std::size_t i = 0; i < buff.size(); ++i) {
    strbuff.insert(i * str_length, buff[i]);
  }
  auto dspace = this->select(offset, buff.size());
  dspace.write_raw(strbuff.data(), dspace.getDataType());

  return offset + buff.size();
}

std::size_t Dataset::write(const internal::VariantBuffer &vbuff, std::size_t offset,
                           bool allow_dataset_resize) {
  std::size_t new_offset{};
  std::visit(
      [&](const auto &buff) { new_offset = this->write(buff, offset, allow_dataset_resize); },
      vbuff.get());

  return new_offset;
}

std::size_t Dataset::read(std::string &buff, std::size_t offset) const {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  if (offset >= this->size()) {
    this->throw_out_of_range_excp(offset);
  }

  auto h5type = this->get_h5type();
  const auto str_length = h5type.getSize();
  buff.resize(str_length);
  this->select(offset, 1).read(buff.data(), h5type);

  const auto i1 = std::min(buff.find('\0'), buff.size());
  buff.resize(i1);

  return offset + 1;
}

std::size_t Dataset::write(std::string buff, std::size_t offset, bool allow_dataset_resize) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  if (offset >= this->size()) {
    if (allow_dataset_resize) {
      this->resize(offset + 1);
    } else {
      this->throw_out_of_range_excp(offset);
    }
  }

  auto dspace = this->select(offset);

  const auto str_length = dspace.getDataType().getSize();
  buff.resize(str_length, '\0');

  dspace.write_raw(buff.data(), dspace.getDataType());

  return offset + 1;
}

std::size_t Dataset::write(const internal::GenericVariant &vbuff, std::size_t offset,
                           bool allow_dataset_resize) {
  std::size_t new_offset{};
  std::visit(
      [&](const auto &buff) { new_offset = this->write(buff, offset, allow_dataset_resize); },
      vbuff);

  return new_offset;
}

internal::VariantBuffer Dataset::read_all(std::size_t offset) const {
  return this->read_all<internal::VariantBuffer>(offset);
}

internal::GenericVariant Dataset::read_last() const {
  return this->read_last<internal::GenericVariant>();
}

internal::GenericVariant Dataset::read(std::size_t offset) const {
  return this->read<internal::GenericVariant>(offset);
}

HighFive::DataSetAccessProps Dataset::generate_default_dset_access_props(
    const std::size_t chunk_size, const std::size_t cache_size) {
  assert(chunk_size != 0);
  assert(cache_size != 0);
  // https://docs.hdfgroup.org/hdf5/v1_12/group___d_a_p_l.html#ga104d00442c31714ee073dee518f661f
  constexpr double w0 = 0.75;  // default as of HDF5 v12.1
  const auto num_chunks = std::max(std::size_t(1), cache_size / chunk_size);
  constexpr auto &prime_number_table = internal::prime_number_table;

  auto it =
      std::lower_bound(prime_number_table.begin(), prime_number_table.end(), 100 * num_chunks);
  const auto num_slots = it != prime_number_table.end() ? *it : prime_number_table.back();

  HighFive::DataSetAccessProps props{};
  props.add(HighFive::Caching(num_slots, cache_size, w0));
  return props;
}

HighFive::DataSetCreateProps Dataset::generate_default_dset_create_props(
    std::uint_fast8_t compression_lvl, const std::size_t chunk_size) {
  assert(chunk_size != 0);
  HighFive::DataSetCreateProps props{};
  props.add(HighFive::Deflate(conditional_static_cast<std::uint32_t>(compression_lvl)));
  props.add(HighFive::Chunking(chunk_size / sizeof(std::int32_t)));
  return props;
}

HighFive::DataSet Dataset::create_fixed_str_dataset(RootGroup &root_grp, std::string_view path,
                                                    std::size_t max_str_length, std::size_t max_dim,
                                                    const HighFive::DataSetAccessProps &aprops,
                                                    const HighFive::DataSetCreateProps &cprops) {
  assert(max_str_length != 0);

  const auto [group_name, dataset_name] = parse_uri(path);
  auto group = root_grp().getGroup(group_name);
  if (group.exist(dataset_name)) {
    throw std::runtime_error(fmt::format(FMT_STRING("Dataset at URI \"{}\" already exists"), path));
  }

  auto dspace = HighFive::DataSpace({0}, {max_dim});

  // Unfortunately we have to drop down to the C api to create this kind of dataset at the moment
  auto dtype_id = H5Tcopy(H5T_C_S1);
  H5Tset_cset(dtype_id, H5T_CSET_ASCII);
  H5Tset_size(dtype_id, max_str_length);
  H5Tset_strpad(dtype_id, H5T_STR_NULLPAD);

  const auto hid = H5Dcreate(group.getId(), dataset_name.c_str(), dtype_id, dspace.getId(),
                             H5P_DEFAULT, cprops.getId(), aprops.getId());

  if (hid == H5I_INVALID_HID) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to create dataset at URI \"{}\""), path));
  }

  std::ignore = H5Dclose(hid);
  return group.getDataSet(dataset_name);
}

std::pair<std::string, std::string> Dataset::parse_uri(std::string_view uri) {
  const auto pos = uri.rfind('/');
  if (pos == std::string_view::npos) {
    return std::make_pair(std::string{"/"}, std::string{uri});
  }

  if (pos + 1 == uri.size()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Invalid dataset URI \"{}\": URI ends with '/'"), uri));
  }

  // clang-format off
  return std::make_pair(std::string{uri.substr(0, pos)},
                        std::string{uri.substr(pos + 1)});
  // clang-format on
}

void Dataset::throw_out_of_range_excp(std::size_t offset) const {
  assert(offset >= this->size());

  if (this->empty()) {
    throw std::out_of_range(fmt::format(
        FMT_STRING("Caught an attempt to access an element of dataset {}, which is empty"),
        this->uri(), offset, this->size()));
  }

  throw std::out_of_range(fmt::format(
      FMT_STRING("Caught an attempt to access an element past the end of dataset {} ({} > {})"),
      this->uri(), offset, this->size()));
}

void Dataset::throw_out_of_range_excp(std::size_t offset, std::size_t n) const {
  assert(offset + n >= this->size());

  if (this->empty()) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("Caught an attempt to access one or more element(s) of dataset {}, "
                               "which is empty ([{}, {}])"),
                    this->uri(), offset, offset + n));
  }

  throw std::out_of_range(
      fmt::format(FMT_STRING("Caught an attempt to access one or more element(s) past the end of "
                             "dataset {} ([{}-{}] >= {})"),
                  this->uri(), offset, offset + n, this->size()));
}

HighFive::DataType Dataset::get_h5type() const {
  auto h5type = this->_dataset.getDataType();
  if (h5type.isFixedLenStr()) {
    return h5type;
  }
  if (h5type.isVariableStr()) {
    return HighFive::create_datatype<std::string>();
  }

  if (h5type.getClass() != HighFive::DataTypeClass::Enum) {
    return h5type;
  }

  // Useful to suppress warnings about treating enum datasets as plain int datasets
  const auto is_unsigned = H5Tget_sign(h5type.getId()) == H5T_SGN_NONE;
  auto create_dtype = [&]([[maybe_unused]] auto tunsigned, [[maybe_unused]] auto tsigned) {
    using T1 = decltype(tunsigned);
    using T2 = decltype(tsigned);
    static_assert(std::is_unsigned_v<T1>);
    static_assert(std::is_signed_v<T2>);
    static_assert(sizeof(T1) == sizeof(T2));
    return is_unsigned ? HighFive::create_datatype<T1>() : HighFive::create_datatype<T2>();
  };

  switch (h5type.getSize()) {
    case 1:
      return create_dtype(std::uint8_t{}, std::int8_t{});
    case 2:
      return create_dtype(std::uint16_t{}, std::int16_t{});
    case 4:
      return create_dtype(std::uint32_t{}, std::int32_t{});
    case 8:
      return create_dtype(std::uint64_t{}, std::int64_t{});
    default:
      unreachable_code();
  }
}

}  // namespace coolerpp
