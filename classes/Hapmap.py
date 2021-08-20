
from collections import defaultdict
from utils.definitions import *
from utils.utils import *
from classes import Call

class Hapmap():
    """
    Class that holds genotype data from hapmap file
    """
    def __init__(self, in_hmp_file):
        self.get_data(in_hmp_file)
        
    def get_data(self, in_hmp_file):
        self.name = in_hmp_file
        smdata = defaultdict()
        msdata = defaultdict()
        with open(in_hmp_file) as fh:
            flag_data = 0
            for line in fh:
                line = line.strip()
                entries = line.split(CSV)
                if entries[0] == "rs#":
                    flag_data = 1
                    samples = entries[11:]
                    continue
                if line == "":
                    continue
                if flag_data == 0:
                    continue
                marker = entries[0]
                if sample in BLANK_SAMPLES:
                    continue
                for i in range(11, len(entries)):
                    sample = samples[i]
                    if sample in BLANK_SAMPLES:
                        continue
                    if len(entries[i]) == 1:
                        entries[i] = IUPAC[entries[i]]
                    if len(entries[i] == 2):
                        entries[i] = to_grid(entries[i])
                    call = Call(entries[i])
                    smdata[sample].put_data(marker, call)
                    msdata[marker].put_data(sample, call)
        self.smdata = smdata
        self.msdata = msdata
