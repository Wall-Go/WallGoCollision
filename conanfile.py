from conan import ConanFile
from conan.tools.build import check_min_cppstd
from conan.errors import ConanInvalidConfiguration

class WallGoCollisionRecipe(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeToolchain", "CMakeDeps"

    def requirements(self):
        self.requires("gsl/2.7.1")
        self.requires("hdf5/1.14.3")
        self.requires("pybind11/2.11.1")
        self.requires("muparser/2.3.4")
        
        # llvm-openmp requires compiler.cppstd>=17
        self.requires("llvm-openmp/18.1.8")
        
        """
        if self.settings.compiler.cppstd:
            try:
                check_min_cppstd(self, 17)
                self.requires("llvm-openmp/18.1.8")
            except ConanInvalidConfiguration as e:
                print("\n\n!! Warning from WallGoCollision !!\n"
                      "Installing OpenMP through Conan requires compiler.cppstd >= 17 in your Conan 'default' profile.\n"
                      "This is non-fatal: the build will attempt to use your system-wide OpenMP installation.\n\n"
                      )
        """
        
        
    def build_requirements(self):
        self.tool_requires("cmake/[>=3.24]")
