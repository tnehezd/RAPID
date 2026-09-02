from setuptools import setup, Extension
import glob
import os
import subprocess

# Összegyűtjük a forrásokat és az include mappákat
source_files = glob.glob("src/**/*.c", recursive=True)

include_directories = ["include"] + [
    os.path.join("include", d) for d in os.listdir("include") 
    if os.path.isdir(os.path.join("include", d))
]

# HDF5 útvonalak lekérdezése Homebrew-ból
extra_include_dirs = list(include_directories)
extra_library_dirs = []

try:
    brew_prefix = subprocess.check_output(["brew", "--prefix", "hdf5"], text=True).strip()
    extra_include_dirs.append(os.path.join(brew_prefix, "include"))
    extra_library_dirs.append(os.path.join(brew_prefix, "lib"))
except Exception:
    for base in ["/opt/homebrew", "/usr/local"]:
        inc = os.path.join(base, "include")
        lib = os.path.join(base, "lib")
        if os.path.exists(inc):
            extra_include_dirs.append(inc)
        if os.path.exists(lib):
            extra_library_dirs.append(lib)

rapidsim_ext = Extension(
    name="rapidsim._core",
    sources=source_files,
    include_dirs=extra_include_dirs,
    library_dirs=extra_library_dirs,
    libraries=["hdf5"],
)

setup(
    ext_modules=[rapidsim_ext],
)