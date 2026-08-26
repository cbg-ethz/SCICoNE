from setuptools import setup, find_packages
import re
import sys

if sys.version_info.major != 3:
    raise RuntimeError("SCICoNE requires Python 3")

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

with open("scicone/__init__.py") as f:
    version = re.search(r'__version__ = "([^"]+)"', f.read()).group(1)

test_requirements = [
    "pytest>=4.4",
    "pytest-runner>=5.0",
]

notebook_requirements = [
    "nbconvert>=5.4.0",
    "nbformat>=4.4.0",
    "jupyter>=1.0.0",
    "ipython>=7.1.1",
]

setup(
    name="scicone",
    version=version,
    description="Single-cell copy number calling and event history reconstruction.",
    url="https://github.com/cbg-ethz/SCICoNE",
    license="GPL-3.0-or-later",
    author=["Pedro Falé Ferreira", "Mustafa Anil Tuncel"],
    author_email=["pedro.ferreira@bsse.ethz.ch", "tuncel.manil@gmail.com"],
    packages=find_packages(),
    install_requires=requirements,
    python_requires=">=3.6",
    package_data={
        '': ['bin/*tests*', 'bin/*inference*', 'bin/*breakpoint_detection*', 'bin/*simulation*', 'data/*']
    },
    include_package_data=True,
    test_suite="tests",
    tests_require=test_requirements,
    extras_require={'notebooks': notebook_requirements,
                    'test': test_requirements},
)
