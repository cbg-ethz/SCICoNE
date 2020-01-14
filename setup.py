from setuptools import setup, find_packages

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="scicone",
    version='1.0dev',
    description='Package for single-cell CNV phylogenetic inference',
    author=['Pedro Falé Ferreira'],
    author_email=['pedro.ferreira@bsse.ethz.ch'],
    packages=find_packages(),
    install_requires=requirements,
)
