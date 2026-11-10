from setuptools import setup, find_packages

setup(
    name="elastic-ase",
    version="0.1.0",
    description="Second-order elastic constants calculation package for ASE (Atomic Simulation Environment)",
    author="Antigravity",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "ase",
        "spglib"
    ],
    entry_points={
        "console_scripts": [
            "elastic-ase=elastic_ase.cli:main"
        ]
    },
    python_requires=">=3.8",
)
