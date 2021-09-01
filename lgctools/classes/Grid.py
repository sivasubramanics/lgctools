
from collections import defaultdict
from utils.definitions import BLANK_SAMPLES, CSV
from classes.Data import *


class Grid():
    """
    Class that holds genotype data from Grid file
    """

    def __init__(self, in_grid_file):
        self.get_data(in_grid_file)

    def get_data(self, in_grid_file):
        self.name = in_grid_file
        smdata = defaultdict()
        msdata = defaultdict()
        with open(in_grid_file) as fh:
            flag_data = 0
            for line in fh:
                line = line.strip()
                entries = line.split(CSV)
                if entries[0] == "DNA \\ Assay":
                    flag_data = 1
                    markers = entries[1:]
                    continue
                if line == "":
                    continue
                if flag_data == 0:
                    continue
                sample = entries[0]
                if sample in BLANK_SAMPLES:
                    continue
                for i in range(1, len(entries)):
                    marker = markers[i-1]
                    if not marker in msdata:
                        msdata[marker] = MS(marker)
                    if not sample in smdata:
                        smdata[sample] = SM(sample)
                    call = Call(entries[i])
                    smdata[sample].put_data(marker, call)
                    msdata[marker].put_data(sample, call)
        self.smdata = smdata
        self.msdata = msdata
