"""Analysis of LAMMPS chunk/avg output
Define a objects to handle analysis of LAMMPS chunk/avg outputs.
"""

# imports
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from builtins import object
import numpy as np
import pandas as pd
import os

class ChunkFrame(object):
    """Storage class for a frame's chunk outputs
    """
    def __init__(self, timestep, nchunk, totalcount):
        self.timestep = timestep
        self.nchunk = nchunk
        self.totalcount = totalcount
        self.chunks_df = None
        return

    def __call__(self):
        return self.chunks_df

    def add_chunks_dict(self, chunks_dict):
        self.chunks_df = pd.DataFrame(chunks_dict)

class ChunkAvg(object):
    """Parse, load, and analyze data from fix chunk/avg outputs.
    """
    def __init__(self, filepath):
        self.fp = os.path.abspath(filepath)
        self.frames = self._parse()

    def _parse(self):
        frames = list()
        # First three lines are comments with info about the outputs.
        column_names = list()
        with open(self.fp, 'r') as f:
            n_line = 0
            comments = True
            firstchunk = True
            nchunk = None
            for n_line, line in enumerate(f):
                if comments:
                    comment = line.split()[1:]
                    if n_line == 2:
                        for item in comment:
                            column_names.append(item)
                        n_columns = len(column_names)
                        comments = False
                elif (firstchunk and not comments):
                    words = line.split()
                    timestep = int(words[0])
                    nchunk = int(words[1])
                    totalcount = int(words[2])
                    chunkframe = ChunkFrame(timestep, nchunk, totalcount)
                    firstchunk = False
                    n_readchunk = 0
                    chunkdict = dict()
                    for cn in column_names:
                        chunkdict[cn] = list()
                else:
                    n_readchunk+=1
                    words = line.split()
                    for j,word in enumerate(words):
                        if j == 0:
                            chunkdict[column_names[j]].append(int(word))
                        else:
                            chunkdict[column_names[j]].append(float(word))
                    if n_readchunk == nchunk:
                        firstchunk = True
                        chunkframe.add_chunks_dict(chunkdict)
                        frames.append(chunkframe)
        return frames

    def block_average(self, block_size, start=0):
        blocks = list()
        nframes = len(self.frames)
        for i in range(start, nframes, block_size):
            bavg = self.frames[i]().loc[:,'Coord1':]
            for j in range(i+1, i+block_size-1):
                bavg = bavg + self.frames[j]().loc[:,'Coord1':]
            bavg = bavg/block_size
            blocks.append(bavg)
        nblocks = len(blocks)
        avg = blocks[0]
        for i in range(1,nblocks):
            avg = avg + blocks[i]
        avg = avg/nblocks
        return avg, blocks

    def chunk_densities(self):
        for i in range(len(self.frames)):
            chunk_df = self.frames[i]()
            cd = self.chunk_density(chunk_df)
            self.frames[i].chunks_df = chunk_df.assign(density=cd)
        return

    def chunk_volumes(self):
        for i in range(len(self.frames)):
            chunk_df = self.frames[i]()
            cd = self.chunk_volume(chunk_df)
            self.frames[i].chunks_df = chunk_df.assign(volume=cd)
        return        


    @staticmethod
    def chunk_density(chunk_df):
        bw = chunk_df['Coord1'][1] - chunk_df['Coord1'][0]
        hbw = bw/2.
        starts = chunk_df['Coord1'][:] - hbw
        ends = chunk_df['Coord1'][:] + hbw
        Vs = (4/3)*np.pi*starts**3
        Ve = (4/3)*np.pi*ends**3
        dV = Ve - Vs
        density = chunk_df['Ncount']/dV
        return density.values

    @staticmethod
    def chunk_volume(chunk_df):
        bw = chunk_df['Coord1'][1] - chunk_df['Coord1'][0]
        hbw = bw/2.
        starts = chunk_df['Coord1'][:] - hbw
        ends = chunk_df['Coord1'][:] + hbw
        Vs = (4/3)*np.pi*starts**3
        Ve = (4/3)*np.pi*ends**3
        dV = Ve - Vs
        return dV.values
