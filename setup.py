import os
import subprocess
from setuptools import setup, Extension


def get_boost_include_dirs():
    """Get Boost include directories for different platforms."""
    include_dirs = []
    
    # Check for Homebrew on macOS
    try:
        result = subprocess.run(
            ["brew", "--prefix", "boost"],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            boost_prefix = result.stdout.strip()
            include_dirs.append(os.path.join(boost_prefix, "include"))
    except FileNotFoundError:
        pass
    
    # Common Linux paths
    if os.path.exists("/usr/include/boost"):
        include_dirs.append("/usr/include")
    if os.path.exists("/usr/local/include/boost"):
        include_dirs.append("/usr/local/include")
    
    return include_dirs


boost_include_dirs = get_boost_include_dirs()

sources = [
    "cigram/cmodel/generate_graph.cc",
    "cigram/cmodel/distributions.cc",
    "cigram/cmodel/sample.cc",
    "cigram/cmodel/cmodel.cc",
]

lfr_sources = [
    "cigram/lfr_model/benchm.cpp",
    "cigram/lfr_model/cast.cpp",
    "cigram/lfr_model/cc.cpp",
    "cigram/lfr_model/combinatorics.cpp",
    "cigram/lfr_model/random.cpp",
    "cigram/lfr_model/lfr_model.cc",
]

cmodule = Extension(
    "cigram.cmodel",
    sources=sources,
    include_dirs=boost_include_dirs,
    extra_compile_args=["-O3", "-std=c++17"],
)

lfrmodule = Extension(
    "cigram.lfr_model",
    sources=lfr_sources,
    include_dirs=boost_include_dirs,
    extra_compile_args=["-Wno-undef", "-O3", "-std=c++17"],
)

setup(ext_modules=[lfrmodule, cmodule])
