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
    "cigram/lfr_model/histograms.cpp",
    "cigram/lfr_model/print.cpp",
    "cigram/lfr_model/random.cpp",
    "cigram/lfr_model/lfr_model.cc",
]

cmodule = Extension("cigram.cmodel", sources=sources, extra_compile_args=["-Ofast"])
lfrmodule = Extension("cigram.lfr_model", sources=lfr_sources, extra_compile_args=["-Ofast"])

setup(
    name="cigram",
    version="1.0.0",
    description="Circular Gaussian Random gRAph Model - a generator for synthetic complex networks",
    zip_safe=False,
    author="James Gilbert",
    author_email="jamie.gilbert@azimov.co.uk",
    license="GPL",
    entry_points={
    },
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    url="http://cigram.ico2s.org",
    include_package_data=True,
    packages=["cigram", "cigram.network_opt"],
    ext_modules=[lfrmodule, cmodule],
)
