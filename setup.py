import pathlib
from setuptools import setup, find_packages

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="spresso",
    version="1.6.0",
    description="Streamlit GUI for Quantum Espresso - Fork of xespresso with integrated GUI.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/vsrsousa/spresso",
    author="VSR Sousa",
    author_email="vsrsousa@users.noreply.github.com",
    maintainer="VSR Sousa",
    maintainer_email="vsrsousa@users.noreply.github.com",
    license="GPL",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "ase>=3.22.0",
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "matplotlib>=3.4.0",
        "pandas>=1.3.0",
        "tqdm>=4.62.0",
        "pyyaml>=5.4.0",
        "jsonschema>=3.2.0",
        "tabulate>=0.8.9",
    ],
    extras_require={
        "gui": ["streamlit>=1.28.0", "plotly>=5.17.0", "py3Dmol>=2.0.0"],
        "ssh": ["paramiko>=2.12.0"],  # Optional SSH functionality
    },
    entry_points={
        "console_scripts": [
            "spresso-gui=gui.__main__:main",
            "xespresso-gui=gui.__main__:main",  # Keep backward compatibility
        ],
    },
    python_requires=">=3.8",
)
