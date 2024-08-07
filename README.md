# LAMMPS Post

![img](assets/LAMMPS-Post_logo_100x100px.png)

![Python version badge](https://img.shields.io/badge/python->=3.9-blue.svg)
[![LAMMPS badge](https://img.shields.io/badge/LAMMPS-23Jun2022-darkblue.svg)](https://www.lammps.org/)
[![license](https://img.shields.io/github/license/blakeaw/lammps-post.svg)](LICENSE)
![version](https://img.shields.io/badge/version-0.1.0-orange.svg)
[![release](https://img.shields.io/github/release-pre/blakeaw/lammps-post.svg)](https://github.com/blakeaw/lammps-post/releases/tag/v0.1.0)


**`lammps-post` facilitates Python-based analysis of thermodynamic and chunk averaging outputs from LAMMPS molecular dynamics simulations.** 

## Table of Contents

 1. [Install](#install)
     1. [Dependencies](#dependencies)
     2. [pip install](#pip-install)
     3. [Manual install](#manual-install)
 2. [License](#license)
 3. [Change Log](#change-log)
 4. [Documentation and Usage](#documentation-and-usage)
     1. [Quick Overview](#quick-overview)
     2. [Example](#example)
     3. [units context manager](#units-context-manager)
     4. [Using units with pysb.macros](#using-units-with-pysbmacros)
     5. [Stochastic Simulation Units](#stochastic-simulation-units)
     6. [unit keyword for Parameter](#unit-keyword-for-parameter)
     7. [Accessing Model units](#accessing-all-the-units-for-a-model)
     8. [Additional Examples](#additional-examples)
     9. [Custom Units](#custom-units)
 5. [Contact](#contact)
 6. [Supporting](#supporting)  
 7. [Other Useful Tools](#other-useful-tools)

------

# Install

| **! Note** |
| :--- |
|  lammps-post is still in version zero development so new versions may not be backwards compatible. |

**lammps-post** installs as the `lmp_post` Python package. It is has been developed with Python 3.9 and LAMMPS version 23Jun2022.

### Dependencies

Note that `lammps-post` has the following core dependencies:
   * [numpy](https://numpy.org/) - developed using version 1.19.2.
   * [pandas](https://pandas.pydata.org/) - developed using version 1.2.1. 


### pip install

You can install `lammps-post` version 0.1.0 with `pip` sourced from the GitHub repo:

##### with git installed:

Fresh install:
```
pip install git+https://github.com/blakeaw/lammps-post@v0.1.0
```
Or to upgrade from an older version:
```
pip install --upgrade git+https://github.com/blakeaw/lammps-post@v0.1.0
```
##### without git installed:

Fresh install:
```
pip install https://github.com/blakeaw/lammps-post/archive/refs/tags/v0.1.0.zip
```
Or to upgrade from an older version:
```
pip install --upgrade https://github.com/blakeaw/lammps-post/archive/refs/tags/v0.1.0.zip
```
### Manual install

First, download the repository. Then from the `lammps-post` folder/directory run
```
pip install .
```

------

# License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

------

# Change Log

See: [CHANGELOG](CHANGELOG.md)

------

# Documentation and Usage

## Quick Overview

`lammps-post` is a Python-based utility for parsing, loading, and analyzing data from output files written during LAMMPS molecular dyanmics simulations. Currently, `lammps-post` defines the following objects:   

 * `ThermoOut('path/to/thermodynamic-output-file')`- reads thermodynamic info controlled by the [thermo](https://docs.lammps.org/thermo.html), [thermo_style](https://docs.lammps.org/thermo_style.html), and [thermo_modify](https://docs.lammps.org/thermo_modify.html) commands from a file (e.g., log.lammps) and parses the data into a Pandas `DataFrame` which can be used for subsequent analyses.
    * If the file contains outputs from multiple sequential runs, a particular run can loaded using the `run_number` option: `ThermoOut('path/to/thermodynamic-output-file', run_number=0)`
 * `ChunkAvg('path/to/ave-chunk-output-file')` - read and parse chunk averaged outputs from files generated by the LAMMPS [fix ave/chunk](https://docs.lammps.org/fix_ave_chunk.html) command. The ChunkAvg class will parse ave/chunk output file and extract the data into a set of Pandas DataFrames.

### Examples

See the following notebooks:
  * [ThermoOut-Example](jupyter-notebooks/ThermoOut-Example.ipynb)
  * [ChunkAvg-Example](jupyter-notebooks/ChunkAvg-Example.ipynb)


------

# Contact

 * **Issues** :bug: : Please open a [GitHub Issue](https://github.com/blakeaw/lammps-post/issues) to report any problems/bugs with the code or its execution, or to make any feature requests.

 * **Discussions** :grey_question: : If you have questions, suggestions, or want to discuss anything else related to the project, feel free to use the [lammps-post Discussions](https://github.com/blakeaw/lammps-post/discussions) board.

------

# Supporting

I'm very happy that you've chosen to use __lammps-post__. If you've found it helpful, here are a few ways you can support its ongoing development:

* **Star** :star: : Show your support by starring the [lammps-post GitHub repository](https://github.com/blakeaw/lammps-post). It helps increase the project's visibility and lets others know it's useful. It also benefits my motivation to continue improving the package!
* **Share** :mega: : Sharing `lammps-post` on your social media, forums, or with your network is another great way to support the project. It helps more people discover `lammps-post`, which in turn motivates me to keep developing!
* **Cite** :books: : Citing or mentioning this software in your work, publications, or projects is another valuable way to support it. It helps spread the word and acknowledges the effort put into its development, which is greatly appreciated!