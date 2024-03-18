import io
import os

from setuptools import setup

import pipemake_utils

with io.open("README.rst", "rt", encoding="utf8") as f:
    readme = f.read()

# Read requirements.txt
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = lib_folder + "/requirements.txt"
requirements = []
if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        requirements = f.read().splitlines()

setup(
    name=pipemake_utils.__name__,
    version=pipemake_utils.__version__,
    project_urls={
    #    "Documentation": pipemake_utils.__docs__,
        "Code": pipemake_utils.__code__,
        "Issue tracker": pipemake_utils.__issue__,
    },
    license=pipemake_utils.__license__,
    url=pipemake_utils.__url__,
    description=pipemake_utils.__summary__,
    long_description_content_type="text/x-rst",
    long_description=readme,
    packages=["pipemake_utils"],
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "process-braker=pipemake_utils.process_braker:main",
            "ped-phenotype-file=pipemake_utils.ped_phenotype_file:main",
            "model-inds = pipemake_utils.model_inds:main",
            "filter-gemma = pipemake_utils.filter_gemma:main",
            "model-pop-files = pipemake_utils.model_pop_files:main",
            "plot-pca = pipemake_utils.plot_pca:main",
        ],
    },
    python_requires=">=3.7",
)
