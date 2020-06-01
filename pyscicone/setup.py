from setuptools import setup, find_packages
import sys

if sys.version_info.major != 3:
    raise RuntimeError("SCICoNE requires Python 3")

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

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
    version="1.0",
    description="Single-cell copy number calling and event history reconstruction.",
    author=["Pedro Falé Ferreira", "Mustafa Anil Tuncel"],
    author_email=["pedro.ferreira@bsse.ethz.ch", "tuncel.manil@gmail.com"],
    packages=find_packages(),
    install_requires=requirements,
    dependency_links=[
        "git+git://github.com/anilbey/PhenoGraph.git@7ef72746688217b2feed2f3f8226da2f304c532c#egg=Phenograph"
    ],
    package_data={
        '': ['bin/*tests*', 'bin/*inference*', 'bin/*breakpoint_detection*', 'bin/*simulation*', 'data/*']
    },
    include_package_data=True,
    test_suite="tests",
    tests_require=test_requirements,
    extras_require={'notebooks': notebook_requirements,
                    'test': test_requirements},
)
