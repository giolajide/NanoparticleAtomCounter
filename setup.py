##generated this setup.py using setup-py-cli
##from https://pypi.org/project/setup-py-cli/

from setuptools import setup, find_packages
import io
import os
from NanoparticleAtomCounter import __version__

here = os.path.abspath(os.path.dirname(__file__))

# read the contents of your README.md
with io.open(os.path.join(here, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="NanoparticleAtomCounter",
    version=__version__,
    author="Gbolagade Olajide, Tibor Szilvasi",
    author_email="giolajide@crimson.ua.edu, tszilvasi@crimson.ua.edu",
    description="Estimates atom counts in monometallic nanoparticles"
    " given radius and contact angle",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/szilvasi-group/NanoparticleAtomCounter#readme",
    packages=find_packages(where="."),
    # specify your dependencies here
    python_requires=">=3.10,<3.14", #pyarrow (for streamlit) ends at 3.13. look again someday
    install_requires=[
        "numpy>=1.20,<2.3",
        "pandas>=1.4",
        "streamlit>=1.50,<2",
        "openpyxl>=3.1",
        "xlrd>=2.0",
    ],
    extras_require={
        
        "test": [
            "pytest>=8,<9",
            "ascii-colors",
        ],
        
        "benchmark": [
            "joblib>=1.5,<2",
            "ase>=3.22,<3.27",
            "ascii-colors",
            "tqdm>=4.66,<5",
            "pytest>=8,<9",
        ],

        "dev": [
            "joblib>=1.5,<2",
            "ase>=3.22,<3.27",
            "ascii-colors",
            "tqdm>=4.66,<5",
            "pytest>=8,<9",
        ],
        
    },
    # enable the CLI: `nanoparticle-atom-count`
    entry_points={
        "console_scripts": [
            "nanoparticle-atom-count="
            "NanoparticleAtomCounter.cli.atom_count:main",  # main script
            
            "atom-count-test="
            "NanoparticleAtomCounter.cli.run_tests:main",  # tests script
            
            "atom-count-benchmark="
            "NanoparticleAtomCounter.cli.benchmark:main", # benchmark script
        ],
    },
    classifiers=[
        ##https://pypi.org/classifiers/
        # Project maturity
        "Development Status :: 4 - Beta",
        # Audience & topic
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        # License
        "License :: OSI Approved :: MIT License",
        # Supported Python versions
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    keywords=["nanoparticle", "tem", "electron microscopy", "active site", "catalyst"],
    license="MIT",
    include_package_data=True,
    project_urls={
        "Bug Tracker": "https://github.com/szilvasi-group/NanoparticleAtomCounter/issues",
        "Documentation": "https://github.com/szilvasi-group/NanoparticleAtomCounter#readme",
    },
)
