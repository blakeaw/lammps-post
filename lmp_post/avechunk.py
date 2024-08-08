"""Objects to handle analysis of LAMMPS fix ave/chunk outputs.
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


class AveChunk(object):
    """Load, parse, and analyze data from a file with fix ave/chunk outputs."""

    _int_cols = list(["Chunk", "OrigID"])

    def __init__(self, filepath: str) -> None:
        """Initialize the AveChunk object with a given output file.

        Args:
            filepath (str): Path and name of the output file from which to parse fix ave/chunk data.
        """

        self.fp = os.path.abspath(filepath)
        self.frames = self._parse()
        return

    def __getitem__(self, index: int) -> pd.DataFrame:
        """Returns the corresponding DataFrame of the ave/chunk ouput at a given index.

        Args:
            index (int): The frame index.

        Raises:
            IndexError: If the index value is not an integer.

        Returns:
            pd.DataFrame: The DataFrame at the given index value.
        """
        if isinstance(index, int):
            return self.frames[index]()
        else:
            msg = f"Tried to index with a non-integer index value: {index}"
            raise IndexError(msg)

    def __len__(self) -> int:
        """Returns the number of frames in the AveChunk."""
        return len(self.frames)

    def _parse(self) -> list[ChunkFrame]:
        """Parses the ave/chunk data file.

        Returns:
            list[ChunkFrame]: list of ChunkFrame objects corresponding to
              each time point, or frame, of the input fix ave/chunk output file.
        """
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

    def average(self, start: int = 0, end: int = -1) -> pd.DataFrame:
        """Returns a DataFrame of averaged values over a given set of frames.

        Args:
            start (int, optional): Starting frame index. Defaults to 0.
            end (int, optional): Final frame index. Defaults to -1.

        Returns:
            pd.DataFrame: Averaged values.
        """
        blocks = list()
        nframes = len(self.frames)
        if end == -1:
            end = nframes
        if end > nframes:
            end = nframes
        avg = self[start].loc[:, "Coord1":]
        count = 1
        for i in range(start + 1, end):
            # print(i)
            avg += self[i].loc[:, "Coord1":]
            count += 1
        avg = avg / count
        return avg

    def block_average(
        self, block_size: int, start: int = 0, end: int = -1
    ) -> tuple[pd.DataFrame, list[pd.DataFrame]]:
        """Computes a block average for the given block size.

        Args:
            block_size (int): Size of each block (i.e., number of frames).
            start (int, optional): Starting frame for the averaging. Defaults to 0.
            end (int, optional): Final frame for the averaging. Defaults to -1.

        Returns:
            tuple[pd.DataFrame, list[pd.DataFrame]]: The first item is a DataFrame
              with the overall average values taken over all of the blocks. The second
              item is a list of DataFrames corresponding each individual block's
              average values.
        """
        blocks = list()
        nframes = len(self.frames)
        if end == -1:
            end = nframes
        if end > nframes:
            end = nframes
        # print((end - start) / block_size)
        for i in range(start, end, block_size):
            # print(i)
            bavg = self[i].loc[:, "Coord1":]
            bsize_count = 1
            for j in range(i + 1, i + block_size):
                if j > end - 1:
                    break
                bavg = bavg + self[j].loc[:, "Coord1":]
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

    def binsphere_chunk_densities(self):
        """Compute the density of radial shells for spherically binned chunks.

        This is only valid for radial chunking using compute chunk/atom with the
        bin/sphere style and assumes Coord1 contains the radii of each shell.

        Adds a new \'density\' column to the DataFrame for each time point/frame.
        """
        for i in range(len(self.frames)):
            chunk_df = self[i]
            cd = self._radial_chunk_density(chunk_df)
            self.frames[i].chunks_df = chunk_df.assign(density=cd)
        return

    def binsphere_chunk_volumes(self):
        """Compute the density of radial shells for spherically binned chunks.

        This is only valid for radial chunking using compute chunk/atom with the
        bin/sphere style and assumes Coord1 contains the radii of each shell.

        Adds a new \'volume\' column to the DataFrame for each time point/frame.
        """
        for i in range(len(self.frames)):
            chunk_df = self[i]
            cv = self._radial_chunk_volume(chunk_df)
            self.frames[i].chunks_df = chunk_df.assign(volume=cv)
        return

    def bin3d_chunk_volumes(self, x_delta: float, y_delta: float, z_delta: float):
        """Compute the volume of 3d binned chunks and add it to each DataFrame.

        This is only valid for 3d chunking using compute chunk/atom with the
        bin/3d style.

        Adds a new \'volume\' column to the DataFrame for each time point/frame.

        Args:
            x_delta (float): Thickness of the bins, delta, along the x-direction.
            y_delta (float): Thickness of the bins, delta, along the y-direction.
            z_delta (float): Thickness of the bins, delta, along the z-direction.
        """
        for i in range(len(self.frames)):
            chunk_df = self[i]
            cd = (x_delta * y_delta * z_delta) * np.ones(len(chunk_df["Ncount"]))
            self.frames[i].chunks_df = chunk_df.assign(volume=cd)
        return

    def bin3d_chunk_pressures(self, cid: str):
        """Compute the pressure of 3d binned chunks and add it to each DataFrame.

        This is only valid for 3d chunking using compute chunk/atom with the
        bin/3d style and requires the stress tensor values s_xx, s_yy, and s_zz from
        compute stress/atom be part of the ave/chunk output data.

        Adds a new \'pressure\' column to the DataFrame for each time point/frame.

        Pressure is computed using 
            p = -(c_id[1] + c_id[2] + c_id[3])/(3*volume)
        and corresponds to the average per atom pressure in that chunk (averaged over all
        atoms in the chunk). 

        The bin3d_chunk_volumes function must be called before this one in order to
        compute the volume of each bin.

        Args:
            cid (str): The compute id of the compute stress/atom command.
        """
        df_0 = self[0]
        if "volume" not in df_0.columns:
            raise ValueError(
                "Must call the bin3d_chunk_volumes function first to compute the volume of each chunk before the pressure can be computed."
            )
        for i in range(len(self.frames)):
            chunk_df = self[i]
            cp = self._bin3d_chunk_pressure(chunk_df, cid)
            self.frames[i].chunks_df = chunk_df.assign(pressure=cp)
        return

    def remove_empty_chunks(self, ncount_threshold: int = 0):
        """Removes empty chunks, or those with atom counts lower than a threshold.

        Args:
            ncount_threshold (int, optional): The threshold value for the 
            atom count in a chunk. Defaults to 0.
        """
        for i in range(len(self.frames)):
            chunk_df = self[i]
            is_gz = chunk_df["Ncount"] > ncount_threshold
            chunk_df_gz = chunk_df[is_gz]
            # self.frames[i].chunks_df = chunk_df_gz.reset_index()
            self.frames[i].chunks_df = chunk_df_gz
        return

    @staticmethod
    def _radial_chunk_density(chunk_df):
        bw = chunk_df["Coord1"].values[1] - chunk_df["Coord1"].values[0]
        hbw = bw / 2.0
        starts = chunk_df["Coord1"].values[:] - hbw
        ends = chunk_df["Coord1"].values[:] + hbw
        Vs = (4 / 3) * np.pi * starts**3
        Ve = (4 / 3) * np.pi * ends**3
        dV = Ve - Vs
        density = chunk_df["Ncount"] / dV
        return density.values

    @staticmethod
    def _radial_chunk_volume(chunk_df):
        bw = chunk_df["Coord1"].values[1] - chunk_df["Coord1"].values[0]
        hbw = bw / 2.0
        starts = chunk_df["Coord1"][:] - hbw
        ends = chunk_df["Coord1"][:] + hbw
        Vs = (4 / 3) * np.pi * starts**3
        Ve = (4 / 3) * np.pi * ends**3
        dV = Ve - Vs
        return dV.values

    @staticmethod
    def _bin3d_chunk_pressure(chunk_df, compid):
        ix = "c_{}[1]".format(compid)
        sxx = chunk_df[ix]
        iy = "c_{}[2]".format(compid)
        syy = chunk_df[iy]
        iz = "c_{}[3]".format(compid)
        szz = chunk_df[iz]
        vol = chunk_df["volume"]
        p = -1.0 * (sxx + syy + szz) / (3.0 * vol)
        return p
