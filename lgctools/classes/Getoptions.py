
import  sys
from collections import defaultdict
from utils.definitions import *

class Getoptions():
    
    def __init__(self, arguments):
        self.task = self.get_task(arguments)
        self.arguments = self.get_arguments(arguments)
        
    def get_task(self, arguments):
        task = []
        for arg in arguments:
            if arg in TASKS:
                task.append(arg)
        if not task:
            # sys.stderr.write(f"No task provided in command line.")
            self.print_usage()
            exit(1)
        if len(task) > 1:
            self.print_usage()
            exit(1)
        return task[0]
        
    def get_arguments(self, arguments):
        out_arguments = defaultdict(dict)
        
        TASK_OPTS_REQ[self.task]
        TASK_OPTS_OPTIONAL[self.task]
        
        return out_arguments
    
    
    def print_usage(self):
        usage = USAGE
        print(f"{usage}")
        