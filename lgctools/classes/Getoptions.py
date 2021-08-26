
import  sys
from utils.definitions import *

class Getoptions():
    
    def __init__(self, arguments):
        self.task = self.get_task(arguments)
        # self.arguments = self.get_arguments(arguments)
        
    def get_task(self, arguments):
        task = []
        for arg in arguments:
            if arg in TASKS:
                task.append()
        if not task:
            # sys.stderr.write(f"No task provided in command line.")
            self.print_usage()
        return task
        
    def print_usage(self):
        usage = USAGE
        print(f"{usage}")
        