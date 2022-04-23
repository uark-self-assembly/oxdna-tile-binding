from genericpath import isdir
import pdb
import os
import sys
import time
import json
import getopt
import copy
from shutil import copyfile
from pyoxdna.analysis import analyze_bonds
from utils import JobLauncher, SIM_HOME, EMAIL_ADDRESS, PYOXDNA_HOME, OXDNA_HOME
from pyoxdna.utils import current_time, read_config, write_config
import subprocess
from datetime import datetime
import importlib.util
import random

""" this runs and analyzes ONE simulation at a time """
TESTING = False

#SIM_STEPS_PER_RUN = 200 if TESTING else 50000

#MAX_TRIES = 100

#INPUT_CONF = "/home/$user$/self-assembly/oxDNA-simulations/newer-oxdna-tile-binding/experiments/center_alone/relaxed.conf"
#INPUT_TOP = "/home/$user$/self-assembly/oxDNA-simulations/newer-oxdna-tile-binding/experiments/center_alone/experiment_center_alone.top"

# TODO get an input file that will work
ANALYSIS_INPUT = f"{PYOXDNA_HOME}/../analysis.input"

### bonds attaching tile to structure
##STARTING_BONDS = [
##    # [ID 1, ID 2] where each ID is the id of a nucleotide in the topology file
##    [944, 498],
##    [945, 497],
##    [946, 496],
##    [947, 495],
##    [948, 494]
##]
##
### this program exits when these bonds are formed
##TARGET_BONDS = [
##    [949, 576],
##    [950, 575],
##    [951, 574],
##    [952, 573],
##    [953, 572],
##
##    [954, 456],
##    [955, 455],
##    [956, 454],
##    [957, 453],
##    [958, 452],
##
##    [959, 639],
##    [960, 638],
##    [961, 637],
##    [962, 636],
##    [963, 635],
##]


#def main(input_conf, input_top, sim_conf, bonds_file, num_steps, num_iterations, out_dir, starting_bonds, target_bonds):
def main(sim_conf, num_iterations, starting_bonds, target_bonds, output_dir, sim_conf_dict, seed=1, remove_iter_files=True, num_steps=10000):
    #global SIM_ID
    #global INPUT_CONF
    #global INPUT_TOP
    #global MAX_TRIES

    #conf_file = input_conf
    #top_file = input_top
    tile_still_bound, tile_is_attached = True, False

    #TODO figure out why this isn't being recognized
    try:
        num_steps = sim_conf_dict['steps']
    except KeyError:
        print(f'Using {num_steps} steps from the command line input')
    counter = 0

    orig_sim_conf_dict = copy.deepcopy(sim_conf_dict)

    # if no seed parameter is provided to code, default to the seed utilized by the prior iteration
    if seed is None:
        seed = orig_sim_conf_dict['seed']

    sim_conf_dict['seed'] = seed


    # make output directory
    os.makedirs(output_dir, exist_ok=True)


    analysis_stats_file = os.path.join(output_dir, 'bond_stats.txt')
    with open(analysis_stats_file, "w") as a_s_file:
        a_s_file.write(f'iter,start,orig,new,target_bound_bonds,orig_detached_bonds\n')

    #TEMP HACK
    #sim_dir = '.'

    # Copy topology file into output file for convenience
    topology_file = sim_conf_dict['topology']
    new_topology_file = os.path.join(output_dir, f"{os.path.basename(topology_file)}")
    abs_topology_file = os.path.abspath(new_topology_file)
    copyfile(topology_file, abs_topology_file)
    orig_sim_conf_dict['topology'] = abs_topology_file

    while tile_still_bound and not tile_is_attached and counter < num_iterations:
        print(f"\n{datetime.now()}: Beginning iteration {counter} of {num_iterations-1} (of {num_steps} steps each)")

        # Run simulation with seed, incremented by the counter
        conf_file, trajectory_file = run_one(sim_conf, output_dir, orig_sim_conf_dict, counter, int(seed)+counter)

        abs_conf_file = os.path.abspath(conf_file)
        abs_trajectory_file = os.path.abspath(trajectory_file)

        print(f'{datetime.now()}: Beginning analysis for iteration {counter} using trajectory file {trajectory_file}')
        #print(f"conf_file: {abs_conf_file}, trajectory_file: {abs_trajectory_file}, topology file: {abs_topology_file}")
        # TODO: can the following use the lastconf file instead of trajectory?
        tile_still_bound, tile_is_attached = analyze_simulation(sim_conf, abs_trajectory_file, abs_topology_file, starting_bonds, target_bonds, analysis_stats_file, counter)

        counter += 1

        if remove_iter_files:
            if counter > 2:
                os.remove(os.path.join(output_dir, f'{counter-2}_trajectory.dat'))
                os.remove(os.path.join(output_dir, f'{counter-2}_simulated.conf'))

    if not tile_still_bound:
        print("THE TILE DETACHED!")

    if tile_is_attached:
        print("TILE IS ATTACHED!")

# TODO add function that updates the timestep of trajectories in the full file once file is complete
# def update_trajectories(traj_file):
    


def same_pair(pair1, pair2):
    return get_sorted_pair(pair1) == get_sorted_pair(pair2)

def get_sorted_pair(pair):
    if pair[0] < pair[1]:
        return [pair[1], pair[0]]
    else:
        return pair

def pair_to_string(pair):
    sorted_pair = get_sorted_pair(pair)
    return f"{sorted_pair[0]}:{sorted_pair[1]}"


def analyze_simulation(conf_file, trajectory_file, top_file, starting_bonds, target_bonds, analysis_stats_file, iteration):
    """ analyze the simulation and return two values: tile_still_bound, tile_is_attached
    
        return tile_still_bound if the tile is bound at the end of the trajectory file
        return tile_is_attached if the tile attaches somewhere in the trajectory file
    """

    bond_data = analyze_bonds(conf_file, trajectory_file, top_file, include_starting_bonds=True, oxDNA_dir=OXDNA_HOME)
    bond_events = bond_data['bond_events'] # list of (time, base1, strand1, action, base2, strand2)

    #print(f'bond_events: {bond_events}')

    target_bond_strings = set([pair_to_string(x) for x in  target_bonds])
    starting_bond_strings = set([pair_to_string(x) for x in  starting_bonds])

    #print(f'starting_bond_strings: {starting_bond_strings}')

    bonds_to_log = target_bond_strings.union(starting_bond_strings)
    all_bonds = set()
    ending_bonds = set()

    # assumes that timestamps are increasing
    for event in bond_events:
        action = event[3]
        bond = pair_to_string([event[1], event[-2]])
        
        if bond not in bonds_to_log:
            continue

        if action == "BINDS":
            all_bonds.add(bond)
            ending_bonds.add(bond)
        elif action == "BREAK":
            ending_bonds.discard(bond)

    num_start = len(starting_bonds)
    num_orig = len(starting_bond_strings.intersection(ending_bonds))
    print(f"Number of original bonds: {num_start} -> {num_orig}")

    num_new = len(target_bond_strings.intersection(all_bonds))
    print(f"Number of new bonds: {num_new}")

    target_bonds_bound = target_bond_strings.intersection(all_bonds)
    orig_bonds_detached = starting_bond_strings - starting_bond_strings.intersection(ending_bonds)

    with open(analysis_stats_file, "a") as a_s_file:
        a_s_file.write(f'{iteration},{num_start},{num_orig},{num_new},{target_bonds_bound},{orig_bonds_detached}\n')

    # all starting bonds are found at end of trajectory
    tile_still_bound = len(starting_bond_strings.intersection(ending_bonds)) > 0

    # some of the target bonds bound
    tile_is_attached = len(target_bond_strings.intersection(all_bonds)) == len(target_bond_strings)

    # print(f'Number of total bonds: {len(target_bond_strings)}')

    return tile_still_bound, tile_is_attached


def run_one(sim_conf, output_dir, orig_sim_conf_dict, counter, seed):
    """ 
        run a topology / conf pair in output_dir and return the resulting conf and topology files

        see oxDNA wiki for more details

        conf -- path to conf file
        top  -- path to top file
    """

    b_split_trajectories = True
    
    #output_sub_dir = os.path.join(output_dir, str(counter))

    # name of conf file which is used to find the output conf file
    #conf_basename = os.path.basename(sim_conf_dict['conf_file']).split(".conf")[0]

    sim_conf_dict = copy.deepcopy(orig_sim_conf_dict)

    #log_file = sim_conf_dict['log_file']
    #print(f"log_file: {log_file}, base: {os.path.basename(log_file)}")
    
    sim_conf_dict['log_file'] = os.path.join(output_dir, f"{counter}_{sim_conf_dict['log_file']}")
    sim_conf_dict['energy_file'] = os.path.join(output_dir, f"{counter}_{sim_conf_dict['energy_file']}")

    if b_split_trajectories:
        # full_trajectory_file will combine all trajectories, while trajectory_file will contain just the current iteration for quicker analysis
        full_trajectory_file = os.path.join(output_dir, f"{os.path.basename(sim_conf_dict['trajectory_file'])}")
        trajectory_file = os.path.join(output_dir, f"{counter}_{sim_conf_dict['trajectory_file']}")
    else:
        trajectory_file = os.path.join(output_dir, f"{os.path.basename(sim_conf_dict['trajectory_file'])}")
    sim_conf_dict['trajectory_file'] = trajectory_file

    last_conf_file = os.path.join(output_dir, f"{counter}_{sim_conf_dict['lastconf_file']}")
    sim_conf_dict['lastconf_file'] = last_conf_file

    if counter > 0:
        sim_conf_dict['conf_file'] = os.path.join(output_dir, f"{counter-1}_{os.path.basename(orig_sim_conf_dict['lastconf_file'])}")



    new_sim_conf = os.path.join(output_dir, f'{counter}_{os.path.basename(sim_conf)}')
    write_config(sim_conf_dict, new_sim_conf)

    command = to_simulation_command(new_sim_conf) + f" --seed={seed}"
    print(command)
 
    # return_code = subprocess.run(command, shell=True, stdout=subprocess.DEVNULL).returncode
    # return_code = 1

    # while return_code:
    return_code = subprocess.run(command.split(),stdout=subprocess.DEVNULL).returncode
        # time.sleep(random.random())
        # return_code = subprocess.check_output(command.split(), stderr=subprocess.STDOUT).decode()
    
 
    if return_code != 0 and return_code != None: # did not exit sucessfully
        raise Exception(return_code)
        # print(return_code)

    # wait for simulation to complete, signaled by lastconf file appearing
    while (not os.path.exists(sim_conf_dict['lastconf_file'])):
        print(f'waiting for last conf...')
        time.sleep(1 if TESTING else 5)

    if b_split_trajectories:
        # combine this trajectory with full trajectory
        with open(full_trajectory_file, "a") as full_traj_file:
            with open(trajectory_file, "r") as traj_file:
                full_traj_file.write(traj_file.read())

    # os.remove(sim_conf)
    # os.remove(traj_file)
##    # wait for simulation to complete, signaled by subdirectory appearing
##    subdirs = get_subdirectories(output_dir, False)
##    while(len(subdirs) == 0):
##        print(f'waiting for subdirs... {subdirs}')
##        time.sleep(1 if TESTING else 5)
##        subdirs = get_subdirectories(output_dir, False)

    #if(len(subdirs) != 1):
    #    raise Exception("too many subdirectories")

    #trajectory_file = os.path.join(subdirs[0], "trajectory.dat")

    #last_conf_file = os.path.join(subdirs[0], f"{conf_basename}.conf")

    # rename last conf file so name doesn't grow each time this function is called
    #os.rename(
        # this is automatically created by run_simulation.py
    #    os.path.join(subdirs[0], f"{conf_basename}-simulate.conf"),
    #    last_conf_file
    #)

    return last_conf_file, trajectory_file

# python run_simulation.py --conf experiments/center_alone/relaxed.conf --top experiments/center_alone/experiment_center_alone.top --md --num_steps 100000

##def dict_to_arg_string(conf_dict):
##    arg_str = ''
##    for k in conf_dict:
##        if k == 'salt_concentration':
##            arg_str = f'{arg_str} --{k}=0.666'
##        else:
##            arg_str = f'{arg_str} --{k}={conf_dict[k]}'
##    return arg_str
        
#def to_simulation_command(conf, top, output_dir, steps=None):
def to_simulation_command(sim_conf):
    #return f"python {PYOXDNA_HOME}/../run_simulation.py --conf {conf} --top {top} --output {output_dir} --md --num_steps {steps}"
    return f'{OXDNA_HOME}/build/bin/oxDNA {sim_conf}'

def has_subdirectory(path):
    return len(get_subdirectories(path)) > 0

def get_subdirectories(path, b_only_linked=True):
    if b_only_linked:
        return [os.path.join(path, x) for x in os.listdir(path) if os.path.islink(os.path.join(path, x))]
    else:
        return [os.path.join(path, x) for x in os.listdir(path)]


def launch_self():
    current_file = os.path.abspath(__file__)

    # directory that holds simulation results
    sim_name = f"auto_analysis_{datetime.now().timestamp()}"
    working_dir = os.path.join(SIM_HOME, sim_name)

    launcher = JobLauncher(
        in_file=job_file,
        command=f'python3 {current_file} --dir {working_dir}'
    )
    
    launcher.launch_job(sim_name=sim_name, working_dir=os.path.join(working_dir, "main"))

def get_bonds_to_analyze(bonds_file):
    starting_bonds = None
    target_bonds   = None

    spec = importlib.util.spec_from_file_location('', bonds_file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    starting_bonds = mod.STARTING_BONDS
    target_bonds = mod.TARGET_BONDS

    return starting_bonds, target_bonds

def display_help():
    print("\nThis script runs and analyzes one oxDNA simulation at a time.")
    print("\nUsage:")
    print("-h, --help\t\t\t\t\tDisplays this help message")
    print("-j <job_file>, --job_file=<job_file>\t\tThe job file to run the simulation with.")
    #print("-c <conf_file>, --conf=<conf_file>\t\tInput configuration file (.conf or .dat)")
    #print("-t <top_file>, --top=<top_file>\t\t\tInput topology file (.top)")
    print("-s <sim_conf_file>, --sim_conf=<sim_conf_file>\tConfiguration file for oxDNA MD simulation")
    print("-b <bonds_file>, --bonds=<bonds_file>\t\tFile ")
    #print("-n <num>, --num_steps=<num>\t\t\tNumber of simulation steps to run in between analyses")
    print("-i <num>, --iterations=<num>\t\t\tMax number of simulation + analysis iterations to run")
    print("-o <out_dir>, --out_dir=<out_dir>\t\tPath to the output directory")
    print("-d <num>, --seed=<num>\t\tSeed to use with oxDNA")
#    print("-g <num>, --gpus=<num>\t\t[Optional] Number of GPUs to use (default is 0, and CPU is used)")
#    print("-p <prefix>, --prefix=<prefix>\t\t[Optional] Prefix to use for names of output files (default is name of config file")
#    print("-d <nu>, --debug=<num>\t\t[Optional] Set level for output of debug messages, 0 (least) to 5 (most)")
    print("")
    exit()

if __name__ == '__main__':
    # The next 3 lines are there to preserve original behavior before adding command line arugments for most options
    args = sys.argv[1:]
    if(len(args) == 0):
        launch_self()
##    elif(args[0] == "--dir"):
##        main(args[1])
##    else:
##        sys.exit()

    # Parse the input parameters
    try:
        # starts at the second element of argv since the first one is the script name
        # extraparams are extra arguments passed after all option/keywords are assigned
        # opts is a list containing the pair "option"/"value"
        opts, extraparams = getopt.getopt(sys.argv[1:], "hj:c:t:s:b:n:i:o:d:",
                            ["help", "job_file=", "conf=", "top=", "sim_conf=", "bonds=", "num_steps=", "iterations=", "out_dir=", "seed="])
        #print 'Opts:',opts
        #print 'Extra parameters:',extraparams
    except:
        print("")
        print("Unrecognized option:" + str(sys.argv[1:]))
        display_help()

    if len(opts) == 0:
        display_help()

    job_file            = None
    input_conf          = None
    input_top           = None
    sim_conf            = None
    bonds_file          = None
    num_steps           = None
    num_iterations      = None
    out_dir             = None
    starting_bonds      = None
    target_bonds        = None
    seed                = None
    
    for o,p in opts:
        if o in ['-h', '--help']:
            display_help()
        elif o in ['-j', '--job_file']:
            job_file = p
        elif o in ['-c', '--conf']:
            input_conf = p
        elif o in ['-t', '--top']:
            input_top = p
        elif o in ['-s', '--sim_conf']:
            sim_conf = p
        elif o in ['-b', '--bonds']:
            bonds_file = p
            starting_bonds, target_bonds = get_bonds_to_analyze(bonds_file)
            if (starting_bonds is None) or (target_bonds is None):
                print("ERROR: Invalid file provided to specify the bonds to be analyzed/watched")
                exit()
        elif o in ['-n', '--num_steps']:
            try:
                num_steps = int(p)
            except:
                print(f"\nERROR: Invalid value '{p}' entered for number of simulation steps (should be an integer)")
                display_help()
        elif o in ['-i', '--iterations']:
            try:
                num_iterations = int(p)
            except:
                print(f"\nERROR: Invalid value '{p}' entered for number of iterations (should be an integer)")
                display_help()
        elif o in ['-o', '--out_dir']:
            out_dir = p

        elif o in ['-d', '--seed']:
            seed = p

##    if input_conf is None:
##        print(f"\nERROR: Invalid arguments. Must specify an input configuration file.")
##        display_help()
##    if input_top is None:
##        print(f"\nERROR: Invalid arguments. Must specify an input topology file.")
##        display_help()
##    if job_file is None:
##        print(f"\nERROR: Invalid arguments. Must specify a job file.")
##        display_help()
    if sim_conf is None:
        print(f"\nERROR: Invalid arguments. Must specify a configuration file for oxDNA MD simulation.")
        display_help()
    if bonds_file is None:
        print(f"\nERROR: Invalid arguments. Must specify a file for the bonds to be analyzed/watched.")
        display_help()
##    if num_steps is None:
##        print(f"\nERROR: Invalid arguments. Must specify the number of simulation steps to run in between analyses.")
##        display_help()
    if num_iterations is None:
        print(f"\nERROR: Invalid arguments. Must specify the max number of simulation + analysis iterations to run.")
        display_help()
    if out_dir is None:
        print(f"\nERROR: Invalid arguments. Must specify the path to the output directory.")
        display_help()
    

    sim_conf_dict = read_config(sim_conf)

    print(f'number of keys in configuration dictionary = {len(sim_conf_dict.keys())}')

    main(sim_conf=sim_conf, num_iterations=num_iterations, starting_bonds=starting_bonds, target_bonds=target_bonds, sim_conf_dict=sim_conf_dict, output_dir=out_dir, seed=seed)
    #main(input_conf=input_conf, input_top=input_top, sim_conf=sim_conf, bonds_file=bonds_file, num_steps=num_steps, num_iterations=num_iterations, out_dir=out_dir)
