from conan import ConanFile
from conan.tools.build import check_min_cppstd
from conan.errors import ConanInvalidConfiguration

import os

class WallGoCollisionRecipe(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeToolchain", "CMakeDeps"

    def requirements(self):
        self.requires("gsl/2.7.1")
        self.requires("hdf5/1.14.3")
        self.requires("pybind11/2.11.1")
        self.requires("muparser/2.3.4")
        
        """Require llvm-openmp recipe if env variable is set.
        This version of OMP doesn't seem to work universally with our lib, so keep the option to use hidden.
        """
        ompEnvVar = os.getenv("WALLGO_USE_CONAN_OMP", "0")
        if ompEnvVar != "0":

            # llvm-openmp requires compiler.cppstd>=17
            if self.settings.compiler.cppstd:
                try:
                    check_min_cppstd(self, 17)
                    self.requires("llvm-openmp/18.1.8")
                except ConanInvalidConfiguration as e:
                    raise RuntimeError("\n\n!! Error from WallGoCollision !!\n"
                        "You have set the 'WALLGO_USE_CONAN_OMP' environment variable which downloads and compiles OpenMP through Conan. "
                        "This option requires compiler.cppstd >= 17 in your Conan 'default' profile. Please modify your profile accordingly, "
                        "or set the environment variable to 0 or undefine it.\n\n"
                        )
        
        
    def build_requirements(self):
        self.tool_requires("cmake/[>=3.18]")
