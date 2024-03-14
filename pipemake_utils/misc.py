import os
import argparse

def confirmFile ():

    # Custom action to confirm file exists
    class customAction(argparse.Action):
        def __call__(self, parser, args, value, option_string = None):
            if not os.path.isfile(value): raise IOError(f'Unable to locate file: {value}')
            setattr(args, self.dest, value)

    return customAction

def confirmDir ():

    # Custom action to confirm directory exists
    class customAction(argparse.Action):
        def __call__(self, parser, args, value, option_string = None):
            if not os.path.isdir(value): raise IOError(f'Unable to locate directory: {value}')
            setattr(args, self.dest, value)

    return customAction