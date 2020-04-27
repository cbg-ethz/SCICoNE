from setuptools import setup, find_packages

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="pyscicone",
    version="1.0dev",
    description="Python wrapper for SCICoNE, a C++ program to infer copy number profiles of single cells and their evolutionary history from WGS.",
    author=["Pedro Falé Ferreira", "Mustafa Anil Tuncel", ],
    author_email=["pedro.ferreira@bsse.ethz.ch", "tuncel.manil@gmail.com"],
    packages=find_packages(),
    install_requires=requirements,
    dependency_links=[
        "git+git://github.com/anilbey/PhenoGraph.git@7ef72746688217b2feed2f3f8226da2f304c532c#egg=Phenograph"
    ],
)
