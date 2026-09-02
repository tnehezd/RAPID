from setuptools import setup
from setuptools.command.build_py import build_py
from wheel.bdist_wheel import bdist_wheel
import subprocess
import os
import shutil

class CustomBuildPy(build_py):
    """Custom build command to compile the C simulation engine before packaging."""
    def run(self):
        print("[SETUP] Compiling C simulation engine from source for this architecture...")
        subprocess.check_call(["make", "clean"])
        subprocess.check_call(["make", "all"])
        
        # Copy the compiled binary inside the rapidsim package directory so it gets bundled
        bin_dir = os.path.join("rapidsim", "bin")
        os.makedirs(bin_dir, exist_ok=True)
        src_bin = os.path.join("bin", "simulation")
        dst_bin = os.path.join(bin_dir, "simulation")
        
        if os.path.exists(src_bin):
            shutil.copy(src_bin, dst_bin)
            print(f"[SETUP] Successfully copied binary to {dst_bin}")
        else:
            raise RuntimeError(f"[SETUP] Error: Binary not found at {src_bin} after make all!")
            
        super().run()

class BinaryBdistWheel(bdist_wheel):
    """Force setuptools to build a platform-specific wheel instead of a pure Python (-any) wheel."""
    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

setup(
    # Metadata is fully handled by pyproject.toml
    packages=["rapidsim"],
    package_data={
        "rapidsim": ["bin/simulation"],  # Ensure the compiled binary is included in the package
    },
    include_package_data=True,
    cmdclass={
        'build_py': CustomBuildPy,
        'bdist_wheel': BinaryBdistWheel,
    },
)