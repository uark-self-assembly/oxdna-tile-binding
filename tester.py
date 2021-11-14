import os
import sys
import time
import json
import random
from datetime import datetime
from pyoxdna import pyoxdna, DNARelaxer
from pyoxdna.analysis import analyze_bonds
from pyoxdna.utils import copy_file, rm, mkdir, copy_file
from utils import JobLauncher, SIM_HOME, OXDNA_HOME
import analyze


def main(output_dir, conf_file, top_file, steps):
    
    copy_file(top_file, 'input.top')

    print('Done')



if __name__ == '__main__':
    args = sys.argv[1:]
    
    if len(args) < 3:
        print("""Razor usage: 
    python run_simulation.py conf_file top_file num_simulations steps
    OR
    cd; ./python3.sh run_simulation.py "conf_file top_file num_simulations steps" <--- note the quotes!
Local usage:
    python run_simulation.py run conf_file top_file steps [output_dir]
""")
        sys.exit()

    if args[0] == 'run': # run main program

        try:
            output_dir = args[4]
            mkdir(output_dir)
        except:
            output_dir = '.'

        main(output_dir, conf_file=args[1], top_file=args[2], steps=int(args[3]))
        
    else: # launch program in pbs script
        launch_self(args)



