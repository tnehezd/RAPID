from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import subprocess
import os
import shutil

class BuildSimBinary(build_py):
    """Custom build command to compile the C simulation engine and embed the binary."""
    def run(self):
        print("[SETUP] Compiling C simulation engine binary...")
        
        # Clean and compile the C project using the existing Makefile
        try:
            subprocess.check_call(["make", "clean"])
            subprocess.check_call(["make", "all"])
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"C compilation failed with error code {e.returncode}")
        
        # Create the internal binary directory inside the rapidsim package
        bin_dir = os.path.join("rapidsim", "bin")
        os.makedirs(bin_dir, exist_ok=True)
        
        # Locate the compiled binary (adjust path if your Makefile outputs elsewhere)
        src_bin = "./bin/simulation"
        dst_bin = os.path.join(bin_dir, "simulation")
        
        if os.path.exists(src_bin):
            shutil.copy(src_bin, dst_bin)
            os.chmod(dst_bin, 0o755)
            print(f"[SETUP] Successfully embedded binary into: {dst_bin}")
        else:
            raise FileNotFoundError(f"Expected compiled binary not found at '{src_bin}'.")
            
        super().run()

setup(
    name="rapidsim",
    version="2.2.2",
    description="Python wrapper for the Rapid Simulation C core engine.",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "rapidsim": ["bin/simulation"],
    },
    cmdclass={
        "build_py": BuildSimBinary,
    },
    entry_points={
        "console_scripts": [
            "rapidsim=rapidsim.wrapper:main",
        ],
    },
    python_requires=">=3.8",
)