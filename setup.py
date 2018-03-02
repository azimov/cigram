from setuptools import setup, Extension, find_packages

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

long_description='''
CiGRAM is a generator for random complex networks with community structure and assortative connections.
CiGRAM can be used as a benchmark for community detection algorithms.

This package also includes LFR benchmark models.
'''

with open('requirements.txt') as f:
    requirements = f.read().splitlines()


setup(
    name="cigram",
    version="0.1.4",
    description="Circular Gaussian Random gRAph Model - a generator for synthetic complex networks",
    long_description=long_description,
    zip_safe=False,
    author="James Gilbert",
    install_requires=requirements,
    author_email="jamie.gilbert@azimov.co.uk",
    license="GPL",
    entry_points={
    },
    ext_modules=[lfrmodule, cmodule],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    url="https://github.com/azimov/cigram",
    include_package_data=True,
    packages=find_packages(),
    project_urls=dict(
        documentation='https://cigram.readthedocs.org',
        visualisation='http://cigram.ico2s.org',
    ),
)
