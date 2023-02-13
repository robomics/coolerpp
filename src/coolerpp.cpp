// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/coolerpp.hpp"

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5Utility.hpp>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "coolerpp/attribute.hpp"
#include "coolerpp/bin_table.hpp"
#include "coolerpp/dataset.hpp"
#include "coolerpp/group.hpp"
#include "coolerpp/internal/numeric_utils.hpp"
#include "coolerpp/internal/prime_number_table.hpp"
#include "coolerpp/internal/type_pretty_printer.hpp"
#include "coolerpp/internal/variant_buff.hpp"
#include "coolerpp/utils.hpp"

namespace coolerpp {

void init_mcool(std::string_view file_path, bool force_overwrite) {
  static constexpr std::array<std::uint64_t, 0> buff{};
  init_mcool(file_path, buff.begin(), buff.end(), force_overwrite);
}

Chromosome find_longest_chromosome(const ChromosomeSet &chroms) {
  if (chroms.empty()) {
    throw std::runtime_error("find_longest_chromosome was called on an empty chromosome map");
  }
  return *std::max_element(chroms.begin(), chroms.end(),
                           [&](const Chromosome &chrom1, const Chromosome &chrom2) {
                             return chrom1.size < chrom2.size;
                           });
}

Chromosome find_chromosome_with_longest_name(const ChromosomeSet &chroms) {
  if (chroms.empty()) {
    throw std::runtime_error(
        "find_chromosome_with_longest_name was called on an empty chromosome set");
  }
  return *std::max_element(chroms.begin(), chroms.end(),
                           [&](const Chromosome &chrom1, const Chromosome &chrom2) {
                             return chrom1.name.size() < chrom2.name.size();
                           });
}

File::File(std::string_view uri, unsigned mode, bool validate)
    : _mode(mode),
      _fp(std::make_unique<HighFive::File>(open_file(uri, _mode, validate))),
      _root_group(open_root_group(*_fp, uri)),
      _groups(open_groups(_root_group)),
      _datasets(open_datasets(_root_group)),
      _attrs(read_standard_attributes(_root_group)),
      _pixel_variant(detect_pixel_type(_root_group)),
      _bins(std::make_unique<BinTable>(
          import_chroms(_datasets.at("chroms/name"), _datasets.at("chroms/length"), false),
          this->bin_size())),
      _index(std::make_unique<Index>(import_indexes(_datasets.at("indexes/chrom_offset"),
                                                    _datasets.at("indexes/bin1_offset"),
                                                    chromosomes(), *_bins, *_attrs.nnz, false))) {
  assert(mode == HighFive::File::ReadOnly || mode == HighFive::File::ReadWrite);
  if (validate) {
    this->validate_bins();
  }
}

File File::open_read_only(std::string_view uri, bool validate) {
  return File(uri, HighFive::File::ReadOnly, validate);
}

File::~File() noexcept {
  try {
    this->finalize();
  } catch (const std::exception &e) {
    fmt::print(stderr, FMT_STRING("{}\n"), e.what());
  } catch (...) {
    fmt::print(stderr,
               FMT_STRING("An unknown error occurred while closing file {}. File is likely "
                          "corrupted or incomplete."),
               this->path());
  }
}

File::operator bool() const noexcept { return !!this->_fp; }

void File::open(std::string_view uri, bool validate) {
  *this = File::open_read_only(uri, validate);
}

void File::close() {
  this->finalize();
  *this = File{};
}

std::string File::uri() const {
  if (this->hdf5_path() == "/") {
    return this->path();
  }
  return fmt::format(FMT_STRING("{}::{}"), this->path(), this->hdf5_path());
}

std::string File::hdf5_path() const { return this->_root_group.hdf5_path(); }
std::string File::path() const {
  if (!*this) {
    return "";
  }
  return this->_fp->getName();
}

std::uint32_t File::bin_size() const noexcept { return this->_attrs.bin_size; }

auto File::chromosomes() const noexcept -> const ChromosomeSet & {
  return this->bins().chromosomes();
}

auto File::bins() const noexcept -> const BinTable & {
  assert(this->_bins);
  return *this->_bins;
}

internal::NumericVariant File::detect_pixel_type(const RootGroup &root_grp, std::string_view path) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  auto dset = root_grp().getDataSet(std::string{path});
  return internal::read_pixel_variant<internal::NumericVariant>(dset);
}

bool File::check_sentinel_attr(const HighFive::Group &grp) {
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

void File::write_sentinel_attr(HighFive::Group grp) {
  assert(!check_sentinel_attr(grp));

  Attribute::write(grp, internal::SENTINEL_ATTR_NAME, internal::SENTINEL_ATTR_VALUE, true);
  grp.getFile().flush();
}

void File::write_sentinel_attr() { File::write_sentinel_attr(this->_root_group()); }

bool File::check_sentinel_attr() { return File::check_sentinel_attr(this->_root_group()); }

HighFive::File File::open_file(std::string_view uri, unsigned int mode, bool validate) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};
  const auto [file_path, root_grp] = parse_cooler_uri(uri);

  const auto new_file = !std::filesystem::exists(file_path);
  HighFive::File f(file_path, mode);
  if (!validate || new_file) {
    return f;
  }

  const auto status = utils::is_cooler(f, root_grp);
  if (!status) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("\"{}\" does not look like a valid Cooler file:\n"
                               "Validation report:\n{}"),
                    uri, status));
  }

  return f;
}

auto File::open_or_create_root_group(HighFive::File &f, std::string_view uri) -> RootGroup {
  if (f.exist(parse_cooler_uri(uri).group_path)) {
    return open_root_group(f, uri);
  }
  return create_root_group(f, uri);
}

auto File::open_root_group(const HighFive::File &f, std::string_view uri) -> RootGroup {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  return {f.getGroup(parse_cooler_uri(uri).group_path)};
}

auto File::create_root_group(HighFive::File &f, std::string_view uri, bool write_sentinel_attr)
    -> RootGroup {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  auto grp = f.createGroup(parse_cooler_uri(uri).group_path);
  if (write_sentinel_attr) {
    Attribute::write(grp, internal::SENTINEL_ATTR_NAME, internal::SENTINEL_ATTR_VALUE);
    f.flush();
  }

  return {grp};
}

auto File::open_groups(const RootGroup &root_grp) -> GroupMap {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
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

auto File::create_groups(RootGroup &root_grp) -> GroupMap {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  GroupMap groups(MANDATORY_GROUP_NAMES.size() + 1);
  groups.emplace(root_grp.hdf5_path(), Group{root_grp, root_grp()});

  std::transform(MANDATORY_GROUP_NAMES.begin(), MANDATORY_GROUP_NAMES.end(),
                 std::inserter(groups, groups.begin()), [&root_grp](const auto group_name) {
                   const auto name = std::string{group_name};
                   auto group_obj = root_grp().createGroup(std::string{group_name});

                   return std::make_pair(name, Group{root_grp, group_obj});
                 });
  return groups;
}

void File::write_standard_attributes(RootGroup &root_grp, const StandardAttributes &attributes,
                                     bool skip_sentinel_attr) {
  assert(attributes.bin_size != 0);
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
  if (attributes.assembly) {
    Attribute::write(root_grp(), "assembly", *attributes.assembly);
  }
  Attribute::write(root_grp(), "bin-size", attributes.bin_size);
  Attribute::write(root_grp(), "bin-type", *attributes.bin_type);
  Attribute::write(root_grp(), "creation-date", *attributes.creation_date);
  Attribute::write(root_grp(), "format", std::string{COOL_MAGIC});
  Attribute::write(root_grp(), "format-url", *attributes.format_url);
  if (!skip_sentinel_attr) {
    static_assert(internal::SENTINEL_ATTR_NAME == "format-version");
    Attribute::write(root_grp(), "format-version", attributes.format_version);
  }
  Attribute::write(root_grp(), "generated-by", *attributes.generated_by);
  Attribute::write(root_grp(), "metadata", *attributes.metadata);
  Attribute::write(root_grp(), "nbins", *attributes.nbins);
  Attribute::write(root_grp(), "nchroms", *attributes.nchroms);
  Attribute::write(root_grp(), "nnz", *attributes.nnz);
  Attribute::write(root_grp(), "storage-mode", *attributes.storage_mode);
  std::visit([&](const auto sum) { Attribute::write(root_grp(), "sum", sum); }, attributes.sum);
}

auto File::open_datasets(const RootGroup &root_grp, std::string_view weight_dataset) -> DatasetMap {
  DatasetMap datasets(MANDATORY_DATASET_NAMES.size() + 1);

  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
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

auto File::read_standard_attributes(const RootGroup &root_grp, bool initialize_missing)
    -> StandardAttributes {
  auto attrs = initialize_missing ? StandardAttributes::init(0) : StandardAttributes::init_empty();
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};

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

  auto read_sum_optional = [&](bool missing_ok) {
    constexpr std::string_view key = "sum";
    auto &buff = attrs.sum;
    if (!Attribute::exists(root_grp(), key) && missing_ok) {
      return false;
    }

    try {
      auto sumv = Attribute::read(root_grp(), key);
      std::visit(
          [&](auto sum) {
            using T = remove_cvref_t<decltype(sum)>;
            if constexpr (std::is_unsigned_v<T>) {
              buff = conditional_static_cast<std::uint64_t>(sum);
              return;
            }
            if constexpr (std::is_signed_v<T>) {
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

  read_sum_optional(missing_ok);

  return attrs;
}

auto File::import_chroms(const Dataset &chrom_names, const Dataset &chrom_sizes, bool missing_ok)
    -> ChromosomeSet {
  try {
    [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
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

[[nodiscard]] static std::vector<std::uint64_t> import_chrom_offsets(const Dataset &dset,
                                                                     std::size_t expected_size) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};
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

Index File::import_indexes(const Dataset &chrom_offset_dset, const Dataset &bin_offset_dset,
                           const ChromosomeSet &chroms, const BinTableLazy &bin_table,
                           std::uint64_t expected_nnz, bool missing_ok) {
  try {
    if (bin_offset_dset.empty()) {
      assert(chrom_offset_dset.empty());
      if (missing_ok) {
        return Index{bin_table};
      }
      throw std::runtime_error("index datasets are empty");
    }

    if (bin_offset_dset.size() != bin_table.size() + 1) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("failed to import offsets from {}: expected {} offsets, found {}"),
                      bin_offset_dset.hdf5_path(), bin_table.size() + 1, bin_offset_dset.size()));
    }

    const auto chrom_offsets = import_chrom_offsets(chrom_offset_dset, chroms.size() + 1);

    Index idx{bin_table, expected_nnz};

    std::size_t bin_id = 0;
    std::for_each(bin_offset_dset.begin<std::uint64_t>(), bin_offset_dset.end<std::uint64_t>(),
                  [&](std::uint64_t offset) {
                    if (bin_id < bin_table.size()) {
                      idx.set_offset_by_bin_id(bin_id++, offset);
                    } else {
                      // Last bin
                      assert(bin_id == bin_table.size());
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

void File::validate_bins() const {
  try {
    assert(this->_attrs.bin_type == "fixed");
    auto nchroms = this->dataset("bins/chrom").size();
    auto nstarts = this->dataset("bins/start").size();
    auto nends = this->dataset("bins/end").size();
    if (nchroms != nstarts || nchroms != nends) {
      throw std::runtime_error(fmt::format(FMT_STRING("Datasets have inconsistent sizes:\n"
                                                      " - \"bins/chrom\": {}\n"
                                                      " - \"bins/start\": {}\n"
                                                      " - \"bins/end\": {}\n"
                                                      "Expected {}"),
                                           nchroms, nstarts, nends, this->bins().size()));
    }

    const auto &nbins = nchroms;
    if (nbins != this->bins().size()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Expected {} bins, found {}"), this->bins().size(), nchroms));
    }

    auto chrom_it = this->dataset("bins/chrom").begin<std::uint32_t>();
    auto start_it = this->dataset("bins/start").begin<std::uint32_t>();
    auto end_it = this->dataset("bins/end").begin<std::uint32_t>();

    auto last_chrom = this->dataset("bins/chrom").end<std::uint32_t>();
    auto last_start = this->dataset("bins/start").end<std::uint32_t>();
    auto last_end = this->dataset("bins/end").end<std::uint32_t>();

    std::size_t i = 0;
    for (const Bin &bin : this->bins()) {
      if (chrom_it == last_chrom || start_it == last_start || end_it == last_end) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Expected {} bins, found {}"), this->bins().size(), i));
      }

      if (this->chromosomes().at(*chrom_it).name != bin.chrom.name || *start_it != bin.start ||
          *end_it != bin.end) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Bin #{}: expected {}:{}-{}, found {}:{}-{}"), i,
                        this->chromosomes().at(*chrom_it).name, *start_it, *end_it, bin.chrom.name,
                        bin.start, bin.end));
      }
      ++chrom_it;
      ++start_it;
      ++end_it;
      ++i;
    }

  } catch (const HighFive::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Bin table at URI {}/{} is invalid or corrupted: {}"), this->uri(),
                    this->group("bins")().getPath(), e.what()));
  }
}

auto File::attributes() const noexcept -> const StandardAttributes & { return this->_attrs; }

const internal::NumericVariant &File::pixel_variant() const noexcept {
  return this->_pixel_variant;
}

bool File::has_signed_pixels() const noexcept {
  // clang-format off
  return this->has_pixel_of_type<std::int8_t>()  ||
         this->has_pixel_of_type<std::int16_t>() ||
         this->has_pixel_of_type<std::int32_t>() ||
         this->has_pixel_of_type<std::int64_t>();
  // clang-format on
}

bool File::has_unsigned_pixels() const noexcept {
  // clang-format off
  return this->has_pixel_of_type<std::uint8_t>()  ||
         this->has_pixel_of_type<std::uint16_t>() ||
         this->has_pixel_of_type<std::uint32_t>() ||
         this->has_pixel_of_type<std::uint64_t>();
  // clang-format on
}

bool File::has_integral_pixels() const noexcept {
  return this->has_signed_pixels() || this->has_unsigned_pixels();
}

bool File::has_float_pixels() const noexcept {
  // clang-format off
  return this->has_pixel_of_type<float>()  ||
         this->has_pixel_of_type<double>() ||
         this->has_pixel_of_type<long double>();
  // clang-format on
}

void File::flush() { this->_fp->flush(); }

auto File::index() noexcept -> Index & {
  assert(this->_index);
  return *this->_index;
}
auto File::index() const noexcept -> const Index & {
  assert(this->_index);
  return *this->_index;
}

void File::write_attributes(bool skip_sentinel_attr) {
  assert(this->_attrs.nbins == this->bins().size());
  assert(this->_attrs.nchroms == this->chromosomes().size());
  assert(this->_attrs.nnz == this->_datasets.at("pixels/count").size());

  File::write_standard_attributes(this->_root_group, this->_attrs, skip_sentinel_attr);
  this->flush();
  if (skip_sentinel_attr) {
    using T [[maybe_unused]] = remove_cvref_t<decltype(internal::SENTINEL_ATTR_VALUE)>;
    assert(Attribute::read<T>(this->_root_group(), internal::SENTINEL_ATTR_NAME) ==
           internal::SENTINEL_ATTR_VALUE);
    Attribute::write(this->_root_group(), "format-version", this->_attrs.format_version, true);
    this->flush();
  }
}

void File::write_chromosomes() {
  assert(this->_datasets.contains("chroms/name"));
  assert(this->_datasets.contains("chroms/length"));
  assert(!this->chromosomes().empty());

  File::write_chromosomes(this->dataset("chroms/name"), this->dataset("chroms/length"),
                          this->chromosomes().begin(), this->chromosomes().end());

  this->_attrs.nchroms = this->chromosomes().size();
}

auto File::group(std::string_view group_name) -> Group & {
  try {
    return this->_groups.at(std::string{group_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Group \"{}\" does not exists!"), group_name));
  }
}
auto File::group(std::string_view group_name) const -> const Group & {
  try {
    return this->_groups.at(std::string{group_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Group \"{}\" does not exists!"), group_name));
  }
}

auto File::dataset(std::string_view dataset_name) -> Dataset & {
  try {
    if (dataset_name.front() == '/') {
      dataset_name = dataset_name.substr(1);
    }
    return this->_datasets.at(std::string{dataset_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Dataset \"{}\" does not exists!"), dataset_name));
  }
}

auto File::dataset(std::string_view dataset_name) const -> const Dataset & {
  try {
    if (dataset_name.front() == '/') {
      dataset_name = dataset_name.substr(1);
    }
    return this->_datasets.at(std::string{dataset_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Dataset \"{}\" does not exists!"), dataset_name));
  }
}

StandardAttributes StandardAttributes::init_empty() noexcept {
  StandardAttributes attrs{};

  attrs.bin_type.reset();
  attrs.creation_date.reset();
  attrs.format_url.reset();
  attrs.generated_by.reset();
  attrs.assembly.reset();
  attrs.nbins.reset();
  attrs.nchroms.reset();
  attrs.metadata.reset();
  attrs.storage_mode.reset();

  return attrs;
}

void File::write_bin_table() {
  File::write_bin_table(this->dataset("bins/chrom"), this->dataset("bins/start"),
                        this->dataset("bins/end"), this->bins());

  this->_attrs.nbins = this->bins().size();
}
void File::write_bin_table(Dataset &chrom_dset, Dataset &start_dset, Dataset &end_dset,
                           const BinTable &bin_table) {
  assert(!bin_table.empty());

  chrom_dset.write(bin_table.begin(), bin_table.end(), 0, true,
                   [&](const Bin &bin) { return bin_table.chromosomes().get_id(bin.chrom); });

  start_dset.write(bin_table.begin(), bin_table.end(), 0, true,
                   [&](const Bin &bin) { return bin.start; });

  end_dset.write(bin_table.begin(), bin_table.end(), 0, true,
                 [&](const Bin &bin) { return bin.end; });

  assert(chrom_dset.size() == bin_table.size());
  assert(start_dset.size() == bin_table.size());
  assert(end_dset.size() == bin_table.size());
}

void File::write_indexes() {
  this->index().finalize(*this->_attrs.nnz);
  File::write_indexes(this->dataset("indexes/chrom_offset"), this->dataset("indexes/bin1_offset"),
                      this->index());
}

void File::write_indexes(Dataset &chrom_offset_dset, Dataset &bin_offset_dset, const Index &idx) {
  chrom_offset_dset.write(idx.compute_chrom_offsets(), 0, true);

  bin_offset_dset.write(idx.begin(), idx.end(), 0, true);

  assert(chrom_offset_dset.size() == idx.num_chromosomes() + 1);
  assert(bin_offset_dset.size() == idx.size() + 1);
}

void File::finalize() {
  if (!_fp || !_finalize) {
    assert(!_bins == !_fp);
    assert(!_index == !_fp);
    return;
  }

  assert(this->_bins);
  assert(this->_index);
  try {
    this->write_chromosomes();
    this->write_bin_table();

    _index->nnz() = *_attrs.nnz;
    this->write_indexes();
    this->write_attributes();

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while closing file {}: {}\n"
                               "File is likely corrupted or incomplete"),
                    this->path(), e.what()));
  }
}

auto File::get_last_bin_written() const -> Bin {
  const auto &dset = this->dataset("pixels/bin1_id");
  if (dset.empty()) {
    return this->bins().bin_id_to_coords(0);
  }
  const auto bin1_id = dset.read_last<std::uint64_t>();
  return this->bins().bin_id_to_coords(bin1_id);
}

}  // namespace coolerpp
