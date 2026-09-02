from setuptools import setup
from setuptools.command.build_py import build_py
import subprocess
import os
import shutil

class CustomBuildPy(build_py):
    def run(self):
        print("[SETUP] Compiling C simulation engine from source...")
        subprocess.check_call(["make", "clean"])
        subprocess.check_call(["make", "all"])
        
        # Safe bundling inside rapidsim package data directory
        data_dir = os.path.join("rapidsim", "data")
        os.makedirs(data_dir, exist_ok=True)
        src_bin = os.path.join("bin", "simulation")
        dst_bin = os.path.join(data_dir, "simulation")
        
        if os.path.exists(src_bin):
            shutil.copy(src_bin, dst_bin)
            print(f"[SETUP] Successfully bundled binary to {dst_bin}")
        else:
            raise RuntimeError(f"[SETUP] Error: Binary not found at {src_bin} after make all!")
            
        super().run()

# Lazy/Safe bdist_wheel hook configuration
cmdclass = {'build_py': CustomBuildPy}
try:
    from wheel.bdist_wheel import bdist_wheel
    class BinaryBdistWheel(bdist_wheel):
        def finalize_options(self):
            super().finalize_options()
            self.root_is_pure = False
    cmdclass['bdist_wheel'] = BinaryBdistWheel
except ImportError:
    pass

setup(
    packages=["rapidsim"],
    package_data={
        "rapidsim": ["data/simulation"],
    },
    include_package_data=True,
    cmdclass=cmdclass,
)