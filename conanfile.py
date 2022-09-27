# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


from conan import ConanFile
from conan.tools.build import check_min_cppstd

required_conan_version = ">=1.51.3"


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

    generators = "cmake", "cmake_find_package", "cmake_find_package_multi"

    def configure(self):
        self.options["highfive"].with_boost = False
        self.options["highfive"].with_eigen = False
        self.options["highfive"].with_opencv = False
        self.options["highfive"].with_xtensor = False

    def validate(self):
        if self.settings.compiler.get_safe("cppstd"):
            check_min_cppstd(self, 17)

    def imports(self):
        self.copy("license*", dst="licenses", folder=True, ignore_case=True)
