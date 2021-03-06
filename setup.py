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

Installation - currently only tested on linux.
Requires C++ and libboost installed and in library paths.

To install visit http://www.boost.org/

There is no reason this shouldn't compile on windows if the required libs are present.
'''

with open('requirements.txt') as f:
    requirements = f.read().splitlines()


setup(
    name="cigram",
    version="0.1.5",
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
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering',
    ],
    platforms="GNU/Linux, Mac OS X >= 10.7, Microsoft Windows >= 7",
)
