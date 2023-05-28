// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp_tools/cli.hpp"

#include <CLI/CLI.hpp>
#include <cassert>
#include <cstdint>
#include <regex>
#include <string>

#include "coolerpp/coolerpp.hpp"

namespace coolerpp::tools {

class CoolerFileValidator : public CLI::Validator {
 public:
  CoolerFileValidator() : Validator("Cooler") {
    func_ = [](std::string& uri) -> std::string {
      if (!coolerpp::utils::is_cooler(uri)) {
        if (coolerpp::utils::is_multires_file(uri)) {
          return "URI points to a .mcool file: " + uri;
        }
        return "Not a valid Cooler: " + uri;
      }
      return "";
    };
  }
};

static bool starts_with(std::string_view s, std::string_view prefix) noexcept {
  return s.find(prefix) == 0;
}

[[nodiscard]] static std::string str_replace_all(std::string s, const std::regex& pattern,
                                                 const std::string& replacement) {
  while (std::regex_search(s, pattern)) {
    s = std::regex_replace(s, pattern, replacement);
  }
  return s;
}

class Formatter : public CLI::Formatter {
  [[nodiscard]] inline std::string make_option_opts(const CLI::Option* opt) const override {
    if (!opt->get_option_text().empty()) {
      return opt->get_option_text();
    }

    auto str_contains = [](const auto s, const auto query) {
      return s.find(query) != decltype(s)::npos;
    };

    std::string out;
    if (opt->get_type_size() != 0) {
      // Format default values so that the help string reads like: --my-option=17.0
      if (!opt->get_default_str().empty()) {
        if (starts_with(opt->get_type_name(), "FLOAT")) {
          auto s = opt->get_default_str();
          if (s.find('.') == std::string::npos) {
            s += ".0";
          }
          out += fmt::format(FMT_STRING("={}"), s);
        } else {
          out += fmt::format(FMT_STRING("={}"), opt->get_default_str());
        }
      }

      // Format param domain using open/closed interval notation
      const std::regex pattern(" - ");
      if (const auto& t = opt->get_type_name(); str_contains(t, " in ")) {
        const auto p1 = t.find("[", t.find(" in "));
        const auto p2 = t.find("]", t.find(" in "));
        if (p1 != std::string::npos && p2 != std::string::npos && p2 > p1) {
          out += " " + str_replace_all(t.substr(p1, p2), pattern, ", ");
        }
      } else if (str_contains(t, "POSITIVE")) {
        out += " (0, inf)";
      } else if (str_contains(t, "NONNEGATIVE") || str_contains(t, "UINT")) {
        out += " [0, inf)";
      }

      if (opt->get_expected_max() == CLI::detail::expected_max_vector_size) {
        out += " ...";
      } else if (opt->get_expected_min() > 1) {
        out += fmt::format(FMT_STRING(" x {}"), opt->get_expected());
      }

      if (opt->get_required()) {
        out += " REQUIRED";
      }
    }
    if (!opt->get_envname().empty()) {
      out += fmt::format(FMT_STRING(" ({}: {})"), get_label("env"), opt->get_envname());
    }
    if (!opt->get_needs().empty()) {
      out += fmt::format(FMT_STRING(" {}:"), get_label("needs"));
      for (const auto* op : opt->get_needs()) {
        out += fmt::format(FMT_STRING(" {}"), op->get_name());
      }
    }
    if (!opt->get_excludes().empty()) {
      out += fmt::format(FMT_STRING(" {}:"), get_label("excludes"));
      for (const auto* op : opt->get_excludes()) {
        out += fmt::format(FMT_STRING(" {}"), op->get_name());
      }
    }

    return out;
  }
};

inline const auto IsValidCoolerFile = CoolerFileValidator();

Cli::Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) { this->make_cli(); }

Cli::subcommand Cli::get_subcommand() const noexcept { return this->_subcommand; }
std::string_view Cli::get_printable_subcommand() const noexcept {
  return Cli::subcommand_to_str(this->get_subcommand());
}

auto Cli::parse_arguments() -> Config {
  this->_cli.name(this->_exec_name);
  this->_cli.parse(this->_argc, this->_argv);

  try {
    if (this->_cli.get_subcommand("dump")->parsed()) {
      this->_subcommand = subcommand::dump;
    } else if (this->_cli.get_subcommand("load")->parsed()) {
      this->_subcommand = subcommand::load;
    } else if (this->_cli.get_subcommand("merge")->parsed()) {
      this->_subcommand = subcommand::merge;
    } else {
      this->_subcommand = subcommand::help;
    }
  } catch (const CLI::ParseError& e) {
    //  This takes care of formatting and printing error messages (if any)
    this->_exit_code = this->_cli.exit(e);
    return this->_config;
  } catch (const std::exception& e) {
    this->_exit_code = 1;
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "An unexpected error has occurred while parsing CLI arguments: {}. If you see this "
            "message, please file an issue on GitHub"),
        e.what()));

  } catch (...) {
    this->_exit_code = 1;
    throw std::runtime_error(
        "An unknown error occurred while parsing CLI arguments! If you see this message, please "
        "file an issue on GitHub");
  }
  this->validate();

  this->_exit_code = 0;
  return this->_config;
}

int Cli::exit(const CLI::ParseError& e) const { return this->_cli.exit(e); }

std::string_view Cli::subcommand_to_str(subcommand s) noexcept {
  switch (s) {
    case dump:
      return "dump";
    case load:
      return "load";
    case merge:
      return "merge";
    default:
      assert(s == help);
      return "--help";
  }
}

void Cli::make_dump_subcommand() {
  auto& sc = *this->_cli.add_subcommand("dump", "Dump Cooler data to stdout.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(this->_config.index() == 0);
                    this->_config = DumpConfig{};
                  });

  this->_config = DumpConfig{};
  auto& c = std::get<DumpConfig>(this->_config);

  // clang-format off
  sc.add_option(
      "cooler-uri",
      c.uri,
      "Path to a Cooler file (URI syntax supported).")
      ->check(IsValidCoolerFile)
      ->required();

  sc.add_option(
      "-t,--table",
      c.table,
      "Name of the table to dump.\n")
      ->check(CLI::IsMember({"chroms", "bins", "pixels"}))
      ->capture_default_str();

  sc.add_option(
      "-r,--range",
      c.range1,
      "Coordinates of the genomic regions to be dumped following UCSC-style notation (chr1:0-1000).")
      ->capture_default_str();

  sc.add_option(
      "--range2",
      c.range2,
      "Coordinates of the genomic regions to be dumped following UCSC-style notation (chr1:0-1000).")
      ->capture_default_str();

  sc.add_option(
      "-b,--balanced",
      c.balanced,
      "Apply balancing weight to data")
      ->default_str("false");

  sc.add_flag(
      "--join,!--no-join",
      c.join,
      "Output pixels in BG2 format.")
      ->capture_default_str();

  sc.add_option(
      "--weight-type",
      c.weight_type,
      "TODO.")
      ->capture_default_str();

  // clang-format on

  this->_config = std::monostate{};
}

void Cli::make_load_subcommand() {
  auto& sc = *this->_cli.add_subcommand("load", "Build .cool files from interactions in BG2/COO.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(this->_config.index() == 0);
                    this->_config = LoadConfig{};
                  });

  this->_config = LoadConfig{};
  auto& c = std::get<LoadConfig>(this->_config);

  // clang-format off
  sc.add_option(
      "chrom-sizes",
      c.path_to_chrom_sizes,
      "Path to .chrom.sizes file.")
      ->check(CLI::ExistingFile)
      ->required();

  sc.add_option(
      "bin-size",
      c.bin_size,
      "Bin size (bp).")
      ->check(CLI::PositiveNumber)
      ->required();

  sc.add_option(
      "output-uri",
      c.uri,
      "Path to output Cooler (URI syntax supported).")
      ->required();

  sc.add_option(
      "-f,--format",
      c.format,
      "Input format.")
      ->check(CLI::IsMember({"bg2", "coo"}))
      ->capture_default_str();

  sc.add_option(
      "--assembly",
      c.assembly,
      "Assembly name.")
      ->capture_default_str();

  sc.add_flag(
      "--count-as-float",
      c.count_as_float,
      "Interactions are floats.")
      ->capture_default_str();

  sc.add_flag(
      "--assume-assume_sorted,!--no-assume-sorted",
      c.assume_sorted,
      "Assume input files are already assume_sorted.")
      ->capture_default_str();
  // clang-format on

  this->_config = std::monostate{};
}

void Cli::make_merge_subcommand() {
  auto& sc = *this->_cli.add_subcommand("merge", "Merge coolers.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(this->_config.index() == 0);
                    this->_config = MergeConfig{};
                  });

  this->_config = MergeConfig{};
  auto& c = std::get<MergeConfig>(this->_config);

  // clang-format off
  sc.add_option(
      "input-coolers",
      c.input_uris,
      "Path to two or more Cooler files to be merged (URI syntax supported).")
      ->check(IsValidCoolerFile)
      ->expected(2, std::numeric_limits<int>::max())
      ->required();

  sc.add_option(
      "-o,--output-cooler",
      c.output_uri,
      "Output Cooler (URI syntax supported).\n"
      "When not specified, merged interactions will be printed to stdout.");

  sc.add_flag(
      "-f,--force",
      c.force,
      "Force overwrite output cooler.")
      ->capture_default_str();

  sc.add_flag(
      "--floating-point,!--integral",
      c.floating_point,
      "Store pixels as floating-point numbers.")
      ->capture_default_str();

  // clang-format on

  this->_config = std::monostate{};
}

void Cli::make_cli() {
  this->_cli.name(this->_exec_name);
  this->_cli.description("Coolerpp tools.");
  this->_cli.set_version_flag("-V,--version", "coolerpp-tools-0.0.1");
  this->_cli.require_subcommand(1);

  this->make_dump_subcommand();
  this->make_load_subcommand();
  this->make_merge_subcommand();
}

void Cli::validate_dump_subcommand() const {
  assert(this->_cli.get_subcommand("dump")->parsed());
  /*
  std::vector<std::string> errors;
  const auto& c = std::get<DumpConfig>(this->_config);

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
   */
}

void Cli::validate_load_subcommand() const {
  assert(this->_cli.get_subcommand("load")->parsed());
  /*
  std::vector<std::string> errors;
  const auto& c = std::get<DumpConfig>(this->_config);

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
  */
}

void Cli::validate_merge_subcommand() const {
  assert(this->_cli.get_subcommand("merge")->parsed());
  /*
  std::vector<std::string> errors;
  const auto& c = std::get<MergeConfig>(this->_config);

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
  */
}

void Cli::validate() const {
  if (this->_cli.get_subcommand("dump")->parsed()) {
    this->validate_dump_subcommand();
  } else if (this->_cli.get_subcommand("load")->parsed()) {
    this->validate_load_subcommand();
  } else if (this->_cli.get_subcommand("merge")->parsed()) {
    this->validate_merge_subcommand();
  } else {
    assert(false);
  }
}

}  // namespace coolerpp::tools
