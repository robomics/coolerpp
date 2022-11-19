# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


from conan import ConanFile
from conan.tools.build import check_min_cppstd

required_conan_version = ">=1.53.0"


class Coolerpp(ConanFile):
    name = "coolerpp"
    homepage = "https://github.com/robomics/coolerpp"
    license = "MIT"
    author = "Roberto Rossini (roberros@uio.no)"
    settings = "os", "compiler", "build_type", "arch"
    requires = ["catch2/3.2.0",
                "fast_float/3.7.0",
                "fmt/9.1.0",
                "hdf5/1.12.2",
                "highfive/2.4.1",
                "tsl-hopscotch-map/2.3.0",
                "tsl-ordered-map/1.0.0",
                "zlib/1.2.13"]

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
