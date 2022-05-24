import io
import os
from setuptools import setup, find_packages

NAME = "charon"
DESCRIPTION = (
    "Charon: A package for neutrino flux generation from WIMP annihilation/decay"
)
MAINTAINER = "Qinrui Liu"
MAINTAINER_EMAIL = "qliu@icecube.wisc.edu"
URL = "https://github.com/IceCubeOpenSource/DMFlux"
LICENSE = "MIT"

here = os.path.abspath(os.path.dirname(__file__))


def read(path, encoding="utf-8"):
    with io.open(path, encoding=encoding) as f:
        content = f.read()
    return content


def get_install_requirements(path):
    content = read(path)
    requirements = [
        req for req in content.split("\n") if req != "" and not req.startswith("#")
    ]
    return requirements


LONG_DESCRIPTION = read(os.path.join(here, "README.md"))

# Want to read in package version number from __version__.py
about = {}
with io.open(os.path.join(here, "charon", "__version__.py"), encoding="utf-8") as f:
    exec(f.read(), about)
    VERSION = about["__version__"]

INSTALL_REQUIRES = get_install_requirements(os.path.join(here, "requirements.txt"))

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url=URL,
    author=MAINTAINER,
    author_email=MAINTAINER_EMAIL,
    license=LICENSE,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    setup_requires=["setuptools>=34.4.0"],
    package_data={"charon": ["data/*.hdf5", "models/*.dat", "xsec/*.h5", "xsec/.dat"]},
)
