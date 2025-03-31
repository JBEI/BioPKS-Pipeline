#!/usr/bin/env python
from setuptools import setup,find_packages

def parse_requirements(filename):
    with open(filename, 'r') as file:
        return file.read().splitlines()

setup(
    name="BioPKS-Pipeline",

    version="1.0",

    description="BioPKS Pipeline is a combined retrobiosynthesis tool that seamlessly integrates RetroTide and DORAnet."
                "Retrotide is a synthesis planning software for designing chimeric type I polyketide synthases (PKSs)"
                "Meanwhile, DORAnet is a synthesis planning software developed specifically for enzymatic synthesis",

    author="Yash Chainani and Tyler Backman",

    author_email="ychainani@lbl.gov",

    url="https://github.com/JBEI/BioPKS-Pipeline",

    packages = find_packages(),

    install_requires=parse_requirements('requirements.txt'),

    package_dir={},

    package_data={},

    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research/Development",
        "Intended Audience :: Scientific Engineering",
        "Intended Audience :: Application",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
    ],
)
