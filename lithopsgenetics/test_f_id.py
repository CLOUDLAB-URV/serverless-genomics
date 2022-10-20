import shutil
import os
import re
import subprocess as sp

def printl(string):
    '''
    add flush=True to print statement to enable real-time printing to log; printl = print to log
    '''
    print(str(var)+"\t"+string, flush=True)

var=1


printl("test")