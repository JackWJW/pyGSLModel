from setuptools import setup, find_packages

with open("./README.md", "r") as f:
    long_description = f.read()

setup(
    name="pyGSLModel",
    version="0.0.10",
    description="A python package for modeling GSL metabolism and performing transcriptomic integration",
    package_dir={"":"src"},
    packages=find_packages(where="src"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JackWJW/pyGSLModel",
    author="Jack Welland",
    license="GNU",
    install_requires=[
        "requests",
        "cobra",
        "io",
        "pyfastcore",
        "mygene",
        "pandas",
        "matplotlib",
        "seaborn",
        "imatpy",
        "numpy"
    ],
    extra_requires={"dev": ["pytest","twine"]}
)