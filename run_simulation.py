""" this runs a vanilla simulation. No fancy optimization methods are used
    see bottom of the file for usage
    see get_configs() and launch_self() to change simulation parameters
"""

import os
import sys
import time
import json
import random
import getopt
from datetime import datetime
from pyoxdna import pyoxdna, DNARelaxer
from pyoxdna.analysis import analyze_bonds
from pyoxdna.utils import copy_file, rm, mkdir, read_config, get_config, chain_configs
from utils import JobLauncher, SIM_HOME, OXDNA_HOME, EMAIL_ADDRESS, print_message
#from debug import *
import analyze

def relax(output_dir, step1, step2, test_config, max_attempts=-1):
    """ relax a configuration, cleanup when done, and return path to relaxed configuration"""
    relaxation_dir =  output_dir #os.path.join(output_dir, 'relaxation')
    r = DNARelaxer(output_dir=relaxation_dir, autoclean=False)
    r.relax(step1, step2, test_config, max_attempts)

    relaxed_conf = os.path.join(relaxation_dir, step2['lastconf_file'])
    # relaxed_conf = os.path.join(output_dir, step2['lastconf_file'])
    # copy_file(os.path.join(relaxation_dir, step2['lastconf_file']), relaxed_conf)
    print_message(3, f'Relaxed configuration: {relaxed_conf}')

    #rm(relaxation_dir)
    return relaxed_conf


def launch_self(output_dir, input_conf, input_top, name_prefix, num_sims, num_steps, relax_steps, b_relax, b_MD, max_relax_attempts, GPUs):
    """ submits this script to be run on Razor 
    """

    walltime='72:00:00'
    # walltime='3:00:00'

    input_conf   = os.path.abspath(input_conf)
    input_top    = os.path.abspath(input_top)
    current_file = os.path.abspath(__file__)

    relax = " --relax" if b_relax else ""
    MD = " --md" if b_MD else ""
    sim_steps = f' --num_steps={num_steps}' if num_steps is not None else ""
    relax_steps = f' --relax_steps={relax_steps}' if relax_steps is not None else ""
    relax_attempts = f' --max_relax={max_relax_attempts}' if max_relax_attempts > -1 else ""
    #debug_level = f' --debug={get_debug_level()}' if get_debug_level() > 0 else ""
    GPU_use = f' --gpus={GPUs}' if GPUs > 0 else ""
    prefix = f' --prefix={name_prefix}' if name_prefix is not None else ""
    
    command = f'python3 {current_file} --local --conf={input_conf} --top={input_top}{relax}{MD}{relax_steps}{sim_steps}{relax_attempts}{prefix}'
    print_message(3, f"Command to launch simulations is '{command}'")

    launcher = JobLauncher(
        email=EMAIL_ADDRESS,
        queue='gpu72',
        nodes=1,
        ppn='6',
        walltime=walltime,
        command=command
    )
    
    if(output_dir is not None and num_sims > 1):
        raise Exception("when using a specific output_dir, num_sims should be 1")
    
    for i in range(num_sims):

        if(output_dir is not None):
            dir_name = [x for x in output_dir.split("/") if len(x) > 0][-1]
            sim_name = f"{dir_name}_{datetime.now().timestamp()}"
            working_dir = output_dir
        else:
            working_dir = os.path.join(SIM_HOME, sim_name)
            sim_name = f'pyoxdna_{datetime.now().timestamp()}'

        print_message(3, f"Simulation number {i}: Working_dir is '{working_dir}'")

        launcher.launch_job(sim_name=sim_name, working_dir=working_dir)
        time.sleep(1)


"""
    conf_file = <string>
        path to the starting configuration
    steps = <int>
        length of the simulation, in time steps

    [lastconf_file = <path>]
        path to the file where the last configuration will be dumped

    [output_prefix = <string>]
        the name of all output files will be preceded by this prefix, defaults
        to an empty string

"""
def main(output_dir, conf_file, top_file, name_prefix, sim_steps, relax_steps, b_relax, b_MD, max_relax_attempts, GPUs):
    
    config_list = []
    config_name_list = []
    
    user_config = dict()

    if b_relax:
        relax1_config = get_config(user_config, 'relax1')
        relax2_config = get_config(user_config, 'relax2')
        test_config = get_config(user_config, 'test')
        
        if not relax_steps is None:
            relax2_config['steps'] = relax_steps
        # relax2_config['backend'] = 'CUDA' if GPUs > 0 else 'CPU'
        
        config_list = [relax1_config, relax2_config, test_config]
        config_name_list = ['relax1', 'relax2', 'relax_test']
    if b_MD:
        # sim_config = get_config(user_config, 'simulate')
        sim_config = get_config({"print_conf_interval": 100, "T":"45C"}, 'simulate')
        
        sim_config['steps'] = sim_steps
        sim_config['backend'] = 'CUDA' if GPUs > 0 else 'CPU'
        sim_config['print_energy_every'] = sim_steps  # don't print unnecessary output
        
        config_list.append(sim_config)
        config_name_list.append('simulate')

    config_list = chain_configs(config_list, config_name_list, name_prefix, top_file, conf_file)

    #copy_file(top_file, 'input.top')
    #user_config = {'conf_file': conf_file, 'topology': top_file}
    for config, config_name in zip(config_list, config_name_list):
        print_message(5, f"Configuration '{config_name}':")
        print_message(5, str(config))

    if b_relax:
        print_message(3, f"Beginning relaxation using output dir '{output_dir}'...")
        relax_out_config = relax(output_dir, relax1_config, relax2_config, test_config, max_relax_attempts)
        print_message(3, f"Relaxation complete and got config_file {relax_out_config}")

    if b_MD:
        print_message(3, f"Beginning molecular dynamics simulation for {sim_steps}...")
        sim = pyoxdna(output_dir=output_dir)
        sim.run_one(sim_config)
        print_message(3, f"Molecular dynamics simulation complete")

    #MJP: skip analysis for now...
    #time.sleep(60*5) # sleep 5 mins to let data finish writing
    #analyze.trim_strands(input_file, sim.config['trajectory_file'], top_file)
    
    
    # submit a job on Razor to analyze the simulation
    # TODO this doesn't work
    # input_file = os.path.join(output_dir, input_file)
    # trajectory_file = os.path.join(output_dir, sim.config['trajectory_file'])
    # sim_name = os.path.basename(output_dir)
    # analyze.launch_self([input_file, trajectory_file, top_file], sim_name=sim_name)

def display_help():
    print("\nThis script initiates and controls simulations of systems in oxDNA, including relaxation.")
    print("and/or molecular dynamics simulation")
    print("\nUsage:")
    print("-h, --help\t\t\t\tDisplays this help message")
    print("-l, --local\t\t\t\t[Optional] Run simulations locally")
    print("-r, --relax\t\t\t\t[Optional] Run relaxation (before simulation if doing both)")
    print("-m, --md\t\t\t\t[Optional] Run molecular dynamics simulation")
    print("-c <conf_file>, --conf=<conf_file>\tInput configuration file (.conf or .dat)")
    print("-t <top_file>, --top=<top_file>\t\tInput topology file (.top)")
    print("-s <num>, --num_sims=<num>\t\tNumber of simulations to run in parallel (only when not local)")
    print("-n <num>, --num_steps=<num>\t\tNumber of simulation steps to run")
    print("-o <out_dir>, --output=<out_dir>\t[Optional] Path to output directory")
    print("-x <num>, --relax_steps=<num>\t\t[Optional] Number of relaxation steps to attempt (default is 2000)")
    print("-a <num>, --max_relax=<num>\t\t[Optional] Maximum number of relaxation attempts to make (default is unbounded)")
    print("-g <num>, --gpus=<num>\t\t[Optional] Number of GPUs to use (default is 0, and CPU is used)")
    print("-p <prefix>, --prefix=<prefix>\t\t[Optional] Prefix to use for names of output files (default is name of config file")
#    print("-d <nu>, --debug=<num>\t\t[Optional] Set level for output of debug messages, 0 (least) to 5 (most)")
    print("")
    exit()

if __name__ == '__main__':
    # Parse the input parameters
    try:
        # starts at the second element of argv since the first one is the script name
        # extraparams are extra arguments passed after all option/keywords are assigned
        # opts is a list containing the pair "option"/"value"
        opts, extraparams = getopt.getopt(sys.argv[1:], "hlrmc:t:s:n:o:x:a:g:p:",
                            ["help", "local", "relax", "md", "conf=", "top=", "num_sims=", "num_steps=", "output=", "relax_steps=", "max_relax=", "gpus=", "prefix="])
        #print 'Opts:',opts
        #print 'Extra parameters:',extraparams
    except:
        print("")
        print("Unrecognized option:" + str(sys.argv[1:]))
        display_help()

    if len(opts) == 0:
        display_help()

    b_local             = False
    b_relax             = False
    b_MD                = False
    input_conf          = None
    input_top           = None
    num_sims            = 1
    num_steps           = None
    output_dir          = None
    num_relax           = None
    max_relax_attempts  = -1
    GPUs                = 0
    prefix              = None
    
    for o,p in opts:
        if o in ['-h', '--help']:
            display_help()
        elif o in ['-l', '--local']:
            b_local = True
        elif o in ['-r', '--relax']:
            b_relax = True
        elif o in ['-m', '--md']:
            b_MD = True
        elif o in ['-c', '--conf']:
            input_conf = p
        elif o in ['-t', '--top']:
            input_top = p
        elif o in ['-s', '--num_sims']:
            try:
                num_sims = int(p)
            except:
                print(f"\nERROR: Invalid value '{p}' entered for number of parallel simulations (should be an integer)")
                display_help()
        elif o in ['-n', '--num_steps']:
            try:
                num_steps = int(p)
            except:
                print(f"\nERROR: Invalid value '{p}' entered for number of simulation steps (should be an integer)")
                display_help()
        elif o in ['-x', '--relax_steps']:
            try:
                num_relax = int(p)
            except:
                print(f"\nERROR: Invalid value '{p}' entered for number of relaxation steps (should be an integer)")
                display_help()
        elif o in ['-a', '--max_relax']:
            try:
                max_relax_attempts = int(p)
            except:
                print(f"\nERROR: Invalid value '{p}' entered for maximum number of relaxation attempts (should be an integer)")
                display_help()
        elif o in ['-o', '--output']:
            output_dir = p
        elif o in ['-g', '--gpus']:
            try:
                GPUs = int(p)
            except:
                print(f"\nERROR: Invalid value '{p}' entered for number of GPUs to use (should be a positive integer)")
        elif o in ['-p', '--prefix']:
            prefix = p
##        elif o in ['d', '--debug']:
##            try:
##                new_debug_level = int(p)
##                assert new_debug_level in range(0,6)
##                set_debug_level(new_debug_level)
##            except:
##                print(f"\nERROR: Invalid value '{p}' entered for debug output level (should be an integer from 0 to 5)")
##                display_help()
        else:
            display_help()

    if input_conf is None:
        print("\nERROR: An input configuration file is required\n")
        display_help()
    elif not os.path.exists(input_conf):
        print(f"\nERROR: The input configuration file '{input_conf}' cannot be found/accessed\n")
        display_help()
    else:
        print(f"Using input configuration file '{input_conf}'")

    if input_top is None:
        print("\nERROR: An input topology file is required\n")
        display_help()
    elif not os.path.exists(input_top):
        print(f"\nERROR: The input topology file '{input_top}' cannot be found/accessed\n")
        display_help()
    else:
        print(f"Using input topology file '{input_top}'")

    if b_relax:
        if not num_relax is None:
            print(f"Relaxation will be performed for a maximum of {num_relax} steps")

        if max_relax_attempts > -1:
            print(f"A maximum of {max_relax_attempts} attempts at relaxation will be made")
        else:
            print("Relaxation attempts will continue until successful")
    if b_MD:
        if num_steps is None:
            print("\nERROR: When performing molecular dynamics simulations, a number of simulation steps must be provided\n")
            display_help()
        print(f"Molecular dynamics simulation(s) will be performed for {num_steps} steps")
    if (not b_relax) and (not b_MD):
        print("\nWARNING: Not performing relaxation or molecular dynamics simulation\n")

    if GPUs == 0:
        print("Simulations will be performed on CPU")
    else:
        print(f"Simulations will be performed on {GPUs} GPU(s)")
        
    if prefix is None or len(prefix) == 0:
        head, tail = os.path.split(input_conf)
        prefix = tail.split('.')[0]
    print(f"Prefix of '{prefix}' will be used for output files")

    if b_local:
        print("Simulations will be run on the local machine")
        if not num_sims == 1:
            print("\nERROR: Only one simulation can be run at a time when running locally")
            display_help()
    else:
        print(f"{num_sims} simulation(s) will be run in parallel")


##MJP: Obsolete format for input parameters
##    args = sys.argv[1:]
##    
##    if len(args) < 3:
##        print("""Razor usage: 
##    python run_simulation.py conf_file top_file num_simulations steps
##    OR
##    cd; ./python3.sh run_simulation.py "conf_file top_file num_simulations steps" <--- note the quotes!
##Local usage:
##    python run_simulation.py run conf_file top_file steps [output_dir]
##""")

    #if args[0] == 'run': # run main program
    if b_local:

        try:
            mkdir(output_dir)
        except:
            output_dir = '.'

        main(output_dir, input_conf, input_top, prefix, num_steps, num_relax, b_relax, b_MD, max_relax_attempts, GPUs)
        
    else: # launch program in pbs script
        launch_self(output_dir, input_conf, input_top, prefix, num_sims, num_steps, num_relax, b_relax, b_MD, max_relax_attempts, GPUs)
