from setuptools import setup, find_packages

with open("./README.md", "r") as f:
    long_description = f.read()

setup(
    name="pyGSLModel",
    version="1.0.0",
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
        "pyfastcore",
        "mygene",
        "pandas",
        "matplotlib",
        "seaborn",
        "imatpy",
        "numpy",
        "huggingface_hub",
        "joblib",
        "tqdm",
        "scikit-learn",
        "skorch",
        "scipy",
        "torch",
        "pyvis"

    ],
    extra_requires={"dev": ["pytest","twine"]},
    python_requires=">=3.10"
)