# Run make in all subdirs and calculate everything that is out of date

import os
import subprocess

mother = os.curdir
dirs = [os.path.join(os.path.abspath(mother), d) for d in os.listdir(mother) if os.path.isdir(d)]
procs = []

for d in dirs:
    os.chdir(d)
    print "\nRUNNING make IN:", os.path.abspath(os.curdir), "\n"
    proc = subprocess.Popen(["make"], shell=True)
    proc.wait()