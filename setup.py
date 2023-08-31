import io
import os

from setuptools import setup

import kp_tools

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
    name=kp_tools.__name__,
    version=kp_tools.__version__,
    project_urls={
        "Documentation": kp_tools.__docs__,
        "Code": kp_tools.__code__,
        "Issue tracker": kp_tools.__issue__,
    },
    license=kp_tools.__license__,
    url=kp_tools.__url__,
    description=kp_tools.__summary__,
    long_description_content_type="text/x-rst",
    long_description=readme,
    packages=["kp_tools"],
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "process-braker=kp_tools.process_braker:main",
        ],
    },
    python_requires=">=3.7",
)
