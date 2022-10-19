import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lammps_output_analysis",
    version="0.1.0",
    python_requires=">=3.9",
    install_requires=[
        "numpy",
        "pandas"
    ],
    author="Blake A. Wilson",
    author_email="blake.wilson@utdallas.edu",
    description="Python tools to parse, post-process, and analyze outputs from LAMMPS commands.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NTBEL/lammps_output_analysis",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
