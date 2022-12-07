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
    """Storage class for a frame's chunk outputs"""

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
    """Parse, load, and analyze data from fix chunk/avg outputs."""

    _int_cols = list(["Chunk", "OrigID"])

    def __init__(self, filepath):
        self.fp = os.path.abspath(filepath)
        self.frames = self._parse()

    def _parse(self):
        frames = list()
        # First three lines are comments with info about the outputs.
        column_names = list()
        with open(self.fp, "r") as f:
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
                elif firstchunk and not comments:
                    words = line.split()
                    timestep = int(words[0])
                    nchunk = int(words[1])
                    totalcount = float(words[2])
                    chunkframe = ChunkFrame(timestep, nchunk, totalcount)
                    firstchunk = False
                    n_readchunk = 0
                    chunkdict = dict()
                    for cn in column_names:
                        chunkdict[cn] = list()
                else:
                    n_readchunk += 1
                    words = line.split()
                    for j, word in enumerate(words):
                        cname = column_names[j]
                        if cname in self._int_cols:
                            chunkdict[column_names[j]].append(int(word))
                        else:
                            chunkdict[column_names[j]].append(float(word))
                    if n_readchunk == nchunk:
                        firstchunk = True
                        chunkframe.add_chunks_dict(chunkdict)
                        frames.append(chunkframe)
        return frames

    def average(self, start=0, end=-1):
        blocks = list()
        nframes = len(self.frames)
        if end == -1:
            end = nframes
        if end > nframes:
            end = nframes
        avg = self.frames[start]().loc[:, "Coord1":]
        count = 1
        for i in range(start + 1, end):
            # print(i)
            avg += self.frames[i]().loc[:, "Coord1":]
            count += 1
        avg = avg / count
        return avg

    def block_average(self, block_size, start=0, end=-1):
        blocks = list()
        nframes = len(self.frames)
        if end == -1:
            end = nframes
        if end > nframes:
            end = nframes
        print((end - start) / block_size)
        for i in range(start, end, block_size):
            # print(i)
            bavg = self.frames[i]().loc[:, "Coord1":]
            bsize_count = 1
            for j in range(i + 1, i + block_size):
                if j > end - 1:
                    break
                bavg = bavg + self.frames[j]().loc[:, "Coord1":]
                bsize_count += 1
            if bsize_count == block_size:
                bavg = bavg / bsize_count
                blocks.append(bavg)
        nblocks = len(blocks)
        avg = blocks[0]
        for i in range(1, nblocks):
            avg = avg + blocks[i]
        avg = avg / nblocks
        for i in range(nblocks):
            blocks[i].dropna(inplace=True)
        return avg, blocks

    def radial_chunk_densities(self):
        for i in range(len(self.frames)):
            chunk_df = self.frames[i]()
            cd = self.radial_chunk_density(chunk_df)
            self.frames[i].chunks_df = chunk_df.assign(density=cd)
        return

    def radial_chunk_volumes(self):
        for i in range(len(self.frames)):
            chunk_df = self.frames[i]()
            cd = self.radial_chunk_volume(chunk_df)
            self.frames[i].chunks_df = chunk_df.assign(volume=cd)
        return

    def bin3d_chunk_volumes(self, x_delta, y_delta, z_delta):
        for i in range(len(self.frames)):
            chunk_df = self.frames[i]()
            cd = (x_delta * y_delta * z_delta) * np.ones(len(chunk_df["Ncount"]))
            self.frames[i].chunks_df = chunk_df.assign(volume=cd)
        return

    def bin3d_chunk_pressures(self, cid):
        for i in range(len(self.frames)):
            chunk_df = self.frames[i]()
            cd = self.bin3d_chunk_pressure(chunk_df, cid)
            self.frames[i].chunks_df = chunk_df.assign(pressure=cd)
        return

    def remove_empty_chunks(self, ncount_threshold=0):
        for i in range(len(self.frames)):
            chunk_df = self.frames[i]()
            is_gz = chunk_df["Ncount"] > ncount_threshold
            chunk_df_gz = chunk_df[is_gz]
            # self.frames[i].chunks_df = chunk_df_gz.reset_index()
            self.frames[i].chunks_df = chunk_df_gz
        return

    @staticmethod
    def radial_chunk_density(chunk_df):
        bw = chunk_df["Coord1"].values[1] - chunk_df["Coord1"].values[0]
        hbw = bw / 2.0
        starts = chunk_df["Coord1"].values[:] - hbw
        ends = chunk_df["Coord1"].values[:] + hbw
        Vs = (4 / 3) * np.pi * starts ** 3
        Ve = (4 / 3) * np.pi * ends ** 3
        dV = Ve - Vs
        density = chunk_df["Ncount"] / dV
        return density.values

    @staticmethod
    def radial_chunk_volume(chunk_df):
        bw = chunk_df["Coord1"].values[1] - chunk_df["Coord1"].values[0]
        hbw = bw / 2.0
        starts = chunk_df["Coord1"][:] - hbw
        ends = chunk_df["Coord1"][:] + hbw
        Vs = (4 / 3) * np.pi * starts ** 3
        Ve = (4 / 3) * np.pi * ends ** 3
        dV = Ve - Vs
        return dV.values

    @staticmethod
    def bin3d_chunk_pressure(chunk_df, compid):
        ix = "c_{}[1]".format(compid)
        sx = chunk_df[ix]
        iy = "c_{}[2]".format(compid)
        sy = chunk_df[iy]
        iz = "c_{}[3]".format(compid)
        sz = chunk_df[iz]
        vol = chunk_df["volume"]
        p = -1.0 * (sx + sy + sz) / (3.0 * vol)
        return p
