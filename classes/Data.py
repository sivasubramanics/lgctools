from collections import defaultdict

class SM():
    def __init__(self, name):
        self.sample_name = name
        self.data = defaultdict()
        
    def put_data(self, marker_name, gt_call):
        self.data[marker_name] = gt_call


class MS():
    def __init__(self, name):
        self.marker_name = name
        self.data = defaultdict()
        
    def put_data(self, sample_name, gt_call):
        self.data[sample_name] = gt_call