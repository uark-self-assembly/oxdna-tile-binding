import os
import sys
import time
import json
from pyoxdna.analysis import analyze_bonds
from utils import JobLauncher, SIM_HOME, EMAIL_ADDRESS, OXDNA_HOME
from pyoxdna.utils import current_time

""" this analyzes a simulation after it has completed """

def trim_strands(input_file, trajectory_file, top_file, output_top='output.top', output_traj='output_trajectory.dat'):
    """ given an input file, trajectory file, and top file, this function creates two new files,
        output_top and output_traj that contain only strands that bind during the simulation.
        This makes it easier to see tile interaction in large simulations
    """    
    bond_data = analyze_bonds(input_file, trajectory_file, top_file, oxDNA_dir=OXDNA_HOME)
    bond_events = bond_data['bond_events'] # list of (time, base1, strand1, action, base2, strand2)

    with open('bond_data.json','w') as f:
        json.dump(bond_events, f)

    if len(bond_events) == 0:
        print('no bonds found in simulation')
        sys.exit()

    strand_nums_to_keep = []
    [strand_nums_to_keep.extend([x[2], x[-1]]) for x in bond_events]
    strand_nums_to_keep = set(strand_nums_to_keep)
    print(strand_nums_to_keep)


    # transfer important strands to new top file
    strands = [] # list of strand lists of (index, base) tuples

    with open(top_file, 'r') as f:
        last_strand = 0
        for i, line in enumerate(f):
            if i == 0: # skip header
                continue

            parts = line.split()
            strand_num = int(parts[0])
            if strand_num in strand_nums_to_keep:
                if strand_num != last_strand:
                    strands.append([])
                    last_strand = strand_num
                
                strands[-1].append([i-1, parts[1]]) # tuple of index, base

    with open(output_top, 'w') as f:
        f.write(f'{len(strands)*len(strands[0])} {len(strands)}\n')
        index = 0
        for i, strand in enumerate(strands):
            for j, (_, base) in enumerate(strand):
                first = index-1 if j != 0 else -1 # attach to previous unless first in strand
                second = index+1 if j != len(strand)-1 else -1 # attach to next unless last in strand
                    
                f.write(f'{i+1} {base} {first} {second}\n')

                index += 1


    # transfer important strands to new traj file
    indices_to_keep = [x[0] for strand in strands for x in strand]
    traj_header = ''
    with open(trajectory_file, 'r') as f:
        with open(output_traj, 'w') as out:
            index = 0
            for line in f:
                if 't = ' in line or 'b = ' in line or 'E = ' in line or line == '\n':
                    out.write(line)
                    index = 0
                else:
                    if index in indices_to_keep:
                        out.write(line)
                    index += 1

def launch_self(args, sim_name=None):
    """ args is a list of arguments:
        [input_file, trajectory_file, top_file]
        OR
        [directory]   
    """

    if len(args) == 4: # python analyze run input trajectory top
        input_file = args[0]
        trajectory_file = args[1]
        top_file = args[2]
        job_file = args[3]

    elif len(args) == 2: # python analyze run directory
        directory = args[0]
        input_file = os.path.join(directory, 'input.params')
        trajectory_file = os.path.join(directory, 'trajectory.dat')
        top_file = os.path.join(directory, 'input.top')
        job_file = args[1]

    else: # invalid usage
        sys.exit()


    walltime='72:00:00'
    # walltime='3:00:00'
    
    current_file = os.path.abspath(__file__)
    sim_name = sim_name or f'analysis_{current_time()}' 

    launcher = JobLauncher(
	in_file=job_file,
        command=f'python3 {current_file} {input_file} {trajectory_file} {top_file}'
    )
    launcher.launch_job(sim_name=sim_name, working_dir=os.path.join(SIM_HOME, sim_name))

if __name__ == '__main__':
    if len(sys.argv) == 4:
        trim_strands(input_file=sys.argv[1], trajectory_file=sys.argv[2], top_file=sys.argv[3])
    else:
        launch_self(sys.argv[1:])
