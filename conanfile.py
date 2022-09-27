# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


from conan import ConanFile
import os
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain, cmake_layout
from conan.tools.build import check_min_cppstd
from conan.tools.files import get, copy

required_conan_version = ">=1.50.0"


class Coolerpp(ConanFile):
    name = "coolerpp"
    homepage = "https://github.com/robomics/coolerpp"
    license = "MIT"
    author = "Roberto Rossini (roberros@uio.no)"
    settings = "os", "compiler", "build_type", "arch"
    requires = ["catch2/3.1.0@#e6317a112d286cae5f349f02d5cda905",
                "fast_float/3.5.1@#63ccdfa6e4dbc05de4bc598258b6a12f",
                "fmt/9.1.0@#78313935d914fe65715ad39687445de6",
                "hdf5/1.12.2@#b01e96ebe1e351ee1d65ae49a347c29c",
                "highfive/2.4.1@#0be5e237f353ec719b2abcc4ada9d9dd",
                "tsl-hopscotch-map/2.3.0@#497d3f41172cefe2df9ac17692c52734",
                "tsl-ordered-map/1.0.0@#dceec5bd0243a5c659939d0295dfee50",
                "zlib/1.2.12@#a30750797caa71bd61bd0a18189caa28"]

    options = {
        "shared": [True, False],
        "fPIC": [True, False]
    }

    default_options = {
        "shared": False,
        "fPIC": True
    }

    def configure(self):
        if self.options.shared:
            del self.options.fPIC

        self.options["highfive"].with_boost = False
        self.options["highfive"].with_eigen = False
        self.options["highfive"].with_opencv = False
        self.options["highfive"].with_xtensor = False

    def layout(self):
        cmake_layout(self, src_folder="src")

    def validate(self):
        if self.settings.compiler.get_safe("cppstd"):
            check_min_cppstd(self, 17)

    def source(self):
        get(self, **self.conan_data["sources"][self.version], destination=self.source_folder,
            strip_root=True)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["ENABLE_TESTING"] = False
        tc.cache_variables["CMAKE_POLICY_DEFAULT_CMP0077"] = "NEW"  # honor BUILD_SHARED_LIBS
        tc.generate()
        tc = CMakeDeps(self)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(self, "LICENSE", dst=os.path.join(self.package_folder, "licenses"),
             src=self.source_folder)
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property("cmake_file_name", "coolerpp")
        self.cpp_info.set_property("cmake_target_name", "coolerpp")
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []
        self.cpp_info.resdirs = []
        self.cpp_info.requires = ["FastFloat::fast_float",
                                  "HighFive::HighFive",
                                  "tsl::hopscotch_map",
                                  "tsl::ordered_map"]

        # TODO: Remove in Conan 2.0
        self.cpp_info.names["cmake_find_package"] = "coolerpp"
        self.cpp_info.names["cmake_find_package_multi"] = "coolerpp"
