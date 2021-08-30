
import sys
from collections import defaultdict
from utils.definitions import *


class Option():
    def __init__(self):
        pass

    def get_option(self):
        return self.option

    def get_input(self):
        return self.input

    def get_default(self):
        return self.default

    def get_help_message(self):
        return self.help_message

    def set_default(self, default):
        self.default = default

    def set_input(self, input):
        self.input = input


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
        print("I'm in..", arguments)
        for i in range(1, len(arguments), 2):
            flag = arguments[i]
            input = arguments[i+1]
            print(flag, input)

        # TASK_OPTS_REQ[self.task]
        # TASK_OPTS_OPTIONAL[self.task]

        return out_arguments

    def print_usage(self):
        usage = USAGE
        print(f"{usage}")
