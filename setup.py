import os
from setuptools import setup, find_packages

# Utility function to read the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup_info = dict(
    name = "SNIRH-WATER-STATION-DISTRIBUTION",
    version = "1.0",
    author = ["Flavio Henrique da Rocha Fonseca", "Avinash Dhiran"],
    author_email = ["fhrochaf@gmail.com", "dhiranavinash@gmail.com"],
    description = ("Tools to retrieve telemtering water station data from HidroWEB-ANA,"
    "\nto retrieve DEM a compute sub-catchments and river network for a target sub-basin,\n"
    "to compute river network as a graph structure, upstream sub-graphs related to water stations,"
    "\n and to analyse water monitoring station distribution."),
    license = "GNU",
    keywords = "telemetering water stations hydrology monitoring Hydroweb SNIRH Brazilian basins river network",
    url = "https://github.com/fhrochaf/SNIRH-Water-Station-Distribution",
    long_description=read('README.md'),
    classifiers=[
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Topic :: Scientific/Engineering :: Hydrology",
        "Topic :: Scientific/Engineering :: GIS"
    ],
    # Package info
    packages=['src'] + ['src.' + pkg for pkg in find_packages('src')],
)

setup(**setup_info)