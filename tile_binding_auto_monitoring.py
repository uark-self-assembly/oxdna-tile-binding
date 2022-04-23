import os
import sys
import time
import subprocess
from math import ceil
from datetime import datetime
from pyoxdna.analysis import analyze_bonds
from utils import JobLauncher, SIM_HOME, EMAIL_ADDRESS

"""
    This script runs an oxDNA simulation for STEPS steps, stopping when at least one of the
    TARGET_BONDS have formed, or all of the STARTING_BONDS have dissolved. Bonds are analyzed
    every ANALYZE_BONDS_INTERVAL steps. INPUT_CONF AND INPUT_TOP are the initial configuration
    and topology files of the oxDNA simulation to run.

    To launch the script, run
        python tile_binding_auto_monitoring.py
"""

INPUT_CONF = "/home/$user$/experiment/newer-oxdna-tile-binding/experiments/center_alone/relaxed.conf"

INPUT_TOP = "/home/$user$/experiment/newer-oxdna-tile-binding/experiments/center_alone/experiment_center_alone.top"

ANALYZE_BONDS_INTERVAL = 10

STEPS = ANALYZE_BONDS_INTERVAL * 100

# list of [ID 1, ID 2] where each ID is the id of a nucleotide in INPUT_TOP
# program exits when all of these bonds have dissolved
STARTING_BONDS = [
    [944, 498],
    [945, 497],
    [946, 496],
    [947, 495],
    [948, 494]
]

# list of [ID 1, ID 2] where each ID is the id of a nucleotide in INPUT_TOP
# program exits when any of these bonds are formed
TARGET_BONDS = [
    [949, 576],
    [950, 575],
    [951, 574],
    [952, 573],
    [953, 572],

    [954, 456],
    [955, 455],
    [956, 454],
    [957, 453],
    [958, 452],

    [959, 639],
    [960, 638],
    [961, 637],
    [962, 636],
    [963, 635],
]




#
# internal variables
#
ANALYSIS_INPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "analysis.input")

# number of seconds between checking if simulation is finished
# scales with steps per sim, with 1 <= x <= 120
MONITOR_SIM_COMPLETION = min(max(ANALYZE_BONDS_INTERVAL / 5000, 1), 120)

def main(sim_dir):
    global INPUT_CONF
    global INPUT_TOP

    max_tries = ceil(STEPS / ANALYZE_BONDS_INTERVAL)

    conf_file = INPUT_CONF
    tile_still_bound, tile_is_attached = True, False
    counter = 0

    while tile_still_bound and not tile_is_attached and counter <= max_tries:
        print(f"\nTIMESTEP {counter * ANALYZE_BONDS_INTERVAL}")

        # for child simulation process
        output_dir = os.path.join(sim_dir, str(counter))
        conf_file, trajectory_file = run_one(conf_file, INPUT_TOP, output_dir)

        tile_still_bound, tile_is_attached = analyze_simulation(trajectory_file, INPUT_TOP)

        counter += 1

    if tile_is_attached:
        print("TILE IS ATTACHED!")
    else:
        print("TILE IS NOT ATTACHED :(")


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


def analyze_simulation(trajectory_file, top_file):
    """ analyze the simulation and return two values: tile_still_bound, tile_is_attached
    
        return tile_still_bound if the tile is bound at the end of the trajectory file
        return tile_is_attached if the tile attaches somewhere in the trajectory file
    """


    bond_data = analyze_bonds(ANALYSIS_INPUT, trajectory_file, top_file, include_starting_bonds=True)
    bond_events = bond_data['bond_events'] # list of (time, base1, strand1, action, base2, strand2)

    target_bond_strings = set([pair_to_string(x) for x in  TARGET_BONDS])
    starting_bond_strings = set([pair_to_string(x) for x in  STARTING_BONDS])

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


    print(f"Number of original bonds: {len(STARTING_BONDS)} -> {len(starting_bond_strings.intersection(ending_bonds))}")

    print(f"Number of new bonds: {len(target_bond_strings.intersection(all_bonds))}")

    # all starting bonds are found at end of trajectory
    tile_still_bound = len(starting_bond_strings.intersection(ending_bonds)) > 0

    # some of the target bonds bound
    tile_is_attached = len(target_bond_strings.intersection(all_bonds)) > 0

    return tile_still_bound, tile_is_attached


def run_one(conf, top, output_dir):
    """ 
        run a topology / conf pair in output_dir and return the resulting conf and topology files

        see oxDNA wiki for more details

        conf -- path to conf file
        top  -- path to top file
    """

    # name of conf file which is used to find the output conf file
    conf_basename = os.path.basename(conf).split(".conf")[0]

    command = to_simulation_command(conf, top, output_dir)

    return_code = subprocess.run(command.split(), stdout=subprocess.DEVNULL).returncode
    
    if return_code != 0: # did not exit sucessfully
        raise Exception(return_code)

    # wait for simulation to complete, signaled by subdirectory appearing
    subdirs = get_subdirectories(output_dir)
    while(len(subdirs) == 0):
        time.sleep(MONITOR_SIM_COMPLETION)
        subdirs = get_subdirectories(output_dir)

    if(len(subdirs) != 1):
        raise Exception("too many subdirectories")

    trajectory_file = os.path.join(subdirs[0], "trajectory.dat")

    last_conf_file = os.path.join(subdirs[0], f"{conf_basename}.conf")

    # rename last conf file so name doesn't grow each time this function is called
    os.rename(
        # this is automatically created by run_simulation.py
        os.path.join(subdirs[0], f"{conf_basename}-simulate.conf"),
        last_conf_file
    )

    return last_conf_file, trajectory_file

# python run_simulation.py --conf experiments/center_alone/relaxed.conf --top experiments/center_alone/experiment_center_alone.top --md --num_steps 100000

def to_simulation_command(conf, top, output_dir, steps=None):
    return f"python /home/$user$/experiment/newer-oxdna-tile-binding/run_simulation.py --conf {conf} --top {top} --output {output_dir} --md --num_steps {steps or ANALYZE_BONDS_INTERVAL}"

def has_subdirectory(path):
    return len(get_subdirectories(path)) > 0

def get_subdirectories(path):
    return [os.path.join(path, x) for x in os.listdir(path) if os.path.islink(os.path.join(path, x))]


def launch_self():
    current_file = os.path.abspath(__file__)

    # directory that holds simulation results
    sim_name = f"auto_analysis_{datetime.now().timestamp()}"
    working_dir = os.path.join(SIM_HOME, sim_name)

    launcher = JobLauncher(
        in_file="example.job",
        command=f'python3 {current_file} --dir {working_dir}'
    )
    
    launcher.launch_job(sim_name=sim_name, working_dir=os.path.join(working_dir, "main"))


if __name__ == '__main__':
    args = sys.argv[1:]

    if(len(args) == 0):
        launch_self()
    elif(args[0] == "--dir"):
        main(args[1])
    else:
        sys.exit()
    

    
