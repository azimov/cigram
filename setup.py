from setuptools import setup, Extension

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

cmodule = Extension("cigram.cmodel", sources=sources, extra_compile_args=["-Ofast"])
lfrmodule = Extension("cigram.lfr_model", sources=lfr_sources, extra_compile_args=["-Wno-undef", "-Ofast"])

setup(
    name="cigram",
    version="0.1.0",
    description="Circular Gaussian Random gRAph Model - a generator for synthetic complex networks",
    zip_safe=False,
    author="James Gilbert",
    author_email="jamie.gilbert@azimov.co.uk",
    license="GPL",
    entry_points={
    },
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    url="https://github.com/azimov/cigram",
    include_package_data=True,
    packages=["cigram", "cigram.network_opt"],
    download_url="https://github.com/azimov/cigram/archive/0.1.0.tar.gz",
)
