from collections import defaultdict
from utils.definitions import *

class SM():
    """
    Class that holds genotype data in sample fast format [sample][marker] = genotype call
    """
    def __init__(self, name):
        self.sample_name = name
        self.data = defaultdict()
        
    def put_data(self, marker_name, gt_call):
        self.data[marker_name] = gt_call

class MS():
    """
    Class that holds genotype data in marker fast format [marker][sample] = genotype call
    """
    def __init__(self, name):
        self.marker_name = name
        self.data = defaultdict()
        
    def put_data(self, sample_name, gt_call):
        self.data[sample_name] = gt_call

class Call():
    """
    Class that holds genotype allele call and its x and y coordinates
    """
    def __init__(self, call):
        self.call = call
        self.x_value = ""
        self.y_value = ""
    
    def get_xvalue(self):
        return self.x_value
    
    def get_yvalue(self):
        return self.y_value
    
    def put_xvalue(self, x_value):
        self.x_value = x_value
        
    def put_yvalue(self, y_value):
        self.y_value = y_value
        
    def __str__(self):
        return str(self.call)
    
    def __repr__(self):
        return str(self.call)
    
    
class Markermetadata():
    """
    Class that marker metadata
    """
    def __init__(self, marker_name):
        self.name = marker_name
        self.allele_x = ""
        self.allele_y = ""
        self.chr = ""
        self.position = ""

    def __repr__(self):
        return self.name

    def put_allele_x(self, allele_x):
        self.allele_x = allele_x

    def get_allele_x(self):
        return self.allele_x + ":" + self.allele_x

    def put_allele_y(self, allele_y):
        self.allele_y = allele_y

    def get_allele_y(self):
        return self.allele_y + ":" + self.allele_y

    def get_allele_xy(self):
        return self.allele_y + ":" + self.allele_x

    def put_alleles(self, lgc_data_line):
        lgc_entries = lgc_data_line.split(CSV)
        if lgc_entries[0] != self.name:
            print(f"Issue with {lgc_entries[0]} marker. Quiting..")
            exit()
        self.put_allele_y(lgc_entries[2])
        self.put_allele_x(lgc_entries[3])

    def get_alleles(self):
        return [self.get_allele_x(), self.get_allele_y(), self.get_allele_xy()]

    def put_chr(self, chr):
        self.chr = chr

    def put_position(self, pos):
        self.position = pos

    def get_chr(self):
        if self.chr:
            return self.chr
        else:
            return 'NA'

    def get_position(self):
        if self.position:
            return self.position
        else:
            return 'NA'

    def put_marker_info(self, marker_info):
        if marker_info:
            self.put_chr(marker_info[0])
            self.put_position(marker_info[1])
    
    
            

