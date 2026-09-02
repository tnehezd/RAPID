from setuptools import setup
from setuptools.command.build_py import build_py
import subprocess
import os
import shutil

class CustomBuildPy(build_py):
    def run(self):
        print("[SETUP] Compiling C simulation engine from source for this architecture...")
        try:
            subprocess.check_call(["make", "clean"])
            subprocess.check_call(["make", "all"])
        except Exception as e:
            raise RuntimeError(f"C compilation failed during package build: {e}")
        
        bin_dir = os.path.join("rapidsim", "bin")
        os.makedirs(bin_dir, exist_ok=True)
        
        src_bin = os.path.join("bin", "simulation")
        dst_bin = os.path.join(bin_dir, "simulation")
        if os.path.exists(src_bin):
            shutil.copy(src_bin, dst_bin)
            os.chmod(dst_bin, 0o755)
        else:
            raise FileNotFoundError(f"Compiled binary not found at {src_bin}")

        super().run()

setup(
    cmdclass={
        "build_py": CustomBuildPy,
    },
)