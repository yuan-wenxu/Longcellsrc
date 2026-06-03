from pathlib import Path
import sysconfig

from setuptools import Extension, setup


ROOT = Path(__file__).parent
PYTHON_INCLUDE = sysconfig.get_paths()["include"]


setup(
    name="longcellsrc",
    version="0.2.0",
    description="Python bindings for Longcellsrc C++ functions",
    package_dir={"": "python"},
    packages=["longcellsrc"],
    ext_modules=[
        Extension(
            name="longcellsrc._core",
            sources=[
                str(ROOT / "python" / "longcellsrc" / "_core.cpp"),
                str(ROOT / "src" / "isoform_core.cpp"),
                str(ROOT / "src" / "longcellsrc_core.cpp"),
                str(ROOT / "src" / "edit.cpp"),
            ],
            include_dirs=[str(ROOT / "src"), PYTHON_INCLUDE],
            language="c++",
            extra_compile_args=["-std=c++17"],
        )
    ],
)
