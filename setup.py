import setuptools
import os

with open("README.rst", "r", encoding="utf-8") as f:
    long_description = f.read()


with open(os.path.join("flux_sensitivity", "version.py")) as f:
    txt = f.read()
    last_line = txt.splitlines()[-1]
    version_string = last_line.split()[-1]
    version = version_string.strip("\"'")


setuptools.setup(
    name="flux_sensitivity_sebastian-achim-mueller",
    version=version,
    description=(
        "Estimate the flux-sensitivity of a point-source either integral "
        "or differential w.r.t. energy."
    ),
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/cherenkov-plenoscope/flux_sensitivity",
    project_urls={
        "Bug Tracker": (
            "https://github.com/cherenkov-plenoscope/"
            "flux_sensitivity/issues"
        ),
    },
    author="Sebastian Achim Mueller",
    author_email="sebastian-achim.mueller@mpi-hd.mpg.de",
    packages=["flux_sensitivity",],
    package_data={
        "flux_sensitivity": [os.path.join("tests", "resources", "*")]
    },
    install_requires=[
        "binning_utils_sebastian-achim-mueller>=0.0.5",
        "propagate_uncertainties_sebastian-achim-mueller>=0.2.2",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
    ],
)
