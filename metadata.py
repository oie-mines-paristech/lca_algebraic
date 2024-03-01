# This is loaded by setuptools to generate python package metadata

import subprocess
from time import time_ns
import re

def run(args) :
    try:
        return subprocess.run(args, stdout=subprocess.PIPE).stdout.decode('utf-8')
    except:
        # When something goes wrong return master
        return "master"

# Get current git branch
curr_branch = run(["git", "branch", "--show-current"])

with open("VERSION", "r") as f:
    VERSION =  re.sub("\s+", "", f.read())
print(VERSION)

if curr_branch != "master" :
    # Setup version name following PEP 440
    VERSION += f".dev{time_ns()}"

