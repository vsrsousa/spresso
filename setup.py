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
    author="Xing Wang (original xespresso), VSR Sousa (spresso fork)",
    author_email="xingwang1991@gmail.com",
    license="GPL",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3.5",
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=["ase", "numpy", "scipy", "matplotlib"],
    extras_require={
        "gui": ["streamlit>=1.28.0", "plotly>=5.17.0", "py3Dmol>=2.0.0"],
    },
    entry_points={
        "console_scripts": [
            "spresso-gui=gui.__main__:main",
            "xespresso-gui=gui.__main__:main",  # Keep backward compatibility
        ],
    },
    python_requires=">=3.5",
)
