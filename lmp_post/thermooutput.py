"""Analysis of LAMMPS thermo output from lammps log or std out capture
Define objects to handle analysis of LAMMPS thermo command outputs.
"""

# imports
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from builtins import object
import numpy as np
import pandas as pd
import os


class ThermoOut(object):
    """Parse, load, and analyze data from fix thermo command outputs."""

    _int_cols = ["Step"]

    def __init__(self, filepath: str, run_number: int = 0):
        """Initialize the ThermoOut object with filepath and run_number to load.

        Args:
            filepath (str): Path and name of the output file (e.g., log.lammps) to be
                be loaded.
            run_number (int, optional): Index of the run for files containing data
                from multiple sequential runs in the same input script. Defaults to 0.
        """
        self.fp = os.path.abspath(filepath)
        self.run_number = run_number
        self.data = self._parse()

    def _parse(self):
        data = dict()
        # Line before start of thermo out is Per MPI ...
        # First line for thermo out are the output quantity names
        # Line after thermo out is Loop time ...
        column_names = list()
        i_run = 0
        with open(self.fp, "r") as f:
            n_line = 0
            pre = True
            nameline = True
            post = False
            for n_line, line in enumerate(f):
                if pre:
                    words = line.split()
                    if len(words) < 2:
                        pass
                    elif (words[0] == "Per") and (words[1] == "MPI"):
                        if i_run == self.run_number:
                            pre = False
                        else:
                            i_run += 1
                    continue
                elif nameline:
                    names = line.split()
                    for item in names:
                        column_names.append(item)
                        n_columns = len(column_names)
                        for cn in column_names:
                            data[cn] = list()
                        nameline = False
                elif not post:
                    words = line.split()
                    if (words[0] == "Loop") and (words[1] == "time"):
                        post = True
                        break
                    for j, word in enumerate(words):
                        cname = column_names[j]
                        if cname in self._int_cols:
                            data[column_names[j]].append(int(word))
                        else:
                            data[column_names[j]].append(float(word))

        return pd.DataFrame(data)
