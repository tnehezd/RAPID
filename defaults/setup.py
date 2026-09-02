from setuptools import setup, Extension

# Itt mondjuk meg, hogy melyik C fájlokat fordítsa le
rapidsim_ext = Extension(
    name="pyrapid._core",  # Ahogy a Pythonból elérhető lesz
    sources=["src/core.c"],  # Ha több .c fájl van, ide sorold fel őket vesszővel elválasztva
    include_dirs=["include"],  # Itt keresi a .h header fájlokat
)

setup(
    ext_modules=[rapidsim_ext],
)