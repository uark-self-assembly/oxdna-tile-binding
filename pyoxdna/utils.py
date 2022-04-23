import os
import sys
import subprocess
import pdb
from math import sqrt
from pprint import pprint
from signal import SIGTERM
from time import sleep, time
from datetime import datetime
from traceback import print_exc
from multiprocessing import Process
from subprocess import run, Popen, CalledProcessError, check_output
from utils.env_vars import PYOXDNA_HOME
from utils.debug import print_message


def combine(*dicts):
    """ combines multiple dictionaries with later dictionaries 
    getting preference (overriding earlier entries) if keys clash """
    output = {}
    for d in dicts:
        output = {**output, **d}
    return output

def rm(filepath):
    try:
        if os.path.isdir(filepath):
            for f in os.listdir(filepath):
                rm(os.path.join(filepath, f))
            os.rmdir(filepath)
        else:
            os.remove(filepath)
    except FileNotFoundError:
        pass
    return


def mkdir(directory):
    """ ensures that a directory exists """
    assert directory is not None, "mkdir directory is None"
    if not os.path.exists(directory):
        os.makedirs(directory)

def write_config(config, path):
    """ writes config (dictionary) to oxDNA input file at path """
    with open(path, 'w') as f:
        print_message(5, f"Writing to config file '{path}':")
        for key in config:
            if config[key] is not None:
                line = key + '=' + str(config[key])
                print_message(5, f'\t{line}')
                f.write(line + '\n')

def read_config(config_file_name):
    """ reads config (dictionary) from input file at path """
    if not os.path.exists(config_file_name):
        pdb.set_trace()
        print(f"Config file '{config_file_name}' cannot be found/accessed")
        return dict()
    else:
        config = dict()
        with open(config_file_name, 'r') as file:
            print_message(5, f"Reading from config file '{config_file_name}':")
            for line in file:
                line = line.strip()
                if len(line) > 0 and not line.startswith('#'):
                    key, value = line.split('=')
                    key = key.strip()
                    value = value.strip()
                    print_message(5, f'\t{key}={value}')
                    config[key] = value
        return config


def chain_configs(config_list, name_list, name_prefix, input_top, input_conf):
    # Arguments:
    #   name_list is a list of stage names, e.g. ['relax1', 'relax2', 'test', 'md-sim']
    #   config_list is a list of config dicts, e.g. those created from 'min.conf', 'relax.conf', 'molecular-dynamics.conf', and 'molecular-dynamics.conf'
    #   name_prefix is the prefix that should begin the names of all output files
    # The first simulation in the list takes input_top and input_conf as inputs.
    # Each stage outputs a last configuration (.conf) and trajectory file (.dat).
    # Each stage takes as input the output configuration of the previous stage.
    # All stages take the same topology file as input.

    assert len(config_list) == len(name_list), f"In chain_conigs, lengths of config_list and name_list don't match: config_list = {config_list}, name_list = {name_list}"
    for name, config in zip(name_list, config_list):
        
        config['topology'] = input_top
        config['conf_file'] = input_conf
        
        out_conf = name_prefix + "-" + name + ".conf"
        out_traj = name_prefix + "-" + name + ".dat"

        config['lastconf_file'] = out_conf
        config['trajectory_file'] = out_traj

        input_conf = out_conf

    return config_list


def append_str(string, obj):
    return f'{string}_{obj}'

def copy_file(from_path, to_path):
    cmd = ['cp', from_path, to_path]
    run(cmd)

def current_time():
    return str(int(round(time() * 1000)))

def kill_pid_file(path):
    with open(path, "r", encoding='utf-8') as f:
        for pid in f.readlines():
            kill(pid)
            

def kill(pid):
    try:
        os.kill(int(pid), SIGTERM)
    except ProcessLookupError:
        pass




def get_config(user_config, config_type):
    """ returns config dictionary for one step of simulation / relaxation 
        config_type is a string
    """
    
    if len(config_type) == 0:
        return user_config

    if config_type == 'relax1':
        min_config = read_config(f'{PYOXDNA_HOME}/configs/min.conf')
        return combine(min_config, user_config)
    elif config_type == 'relax2':
        relax_config = read_config(f'{PYOXDNA_HOME}/configs/relax.conf')
        return combine(relax_config, user_config)
    elif config_type == 'test':
        test_config = read_config(f'{PYOXDNA_HOME}/configs/molecular-dynamics.conf')
        test_config['steps'] = 0
        return combine(test_config, user_config)
    elif config_type == 'simulate':
        md_config = read_config(f'{PYOXDNA_HOME}/configs/molecular-dynamics.conf')
        return combine(md_config, user_config)
    else:
        assert False, f"invalid configuration '{config_type}' requested in get_config"

##    # this must be in order of least importance for combine() to work correctly
##    defaults = [DEFAULT]
##
##    if 'relax' in config_type:
##        defaults.append(RELAX)
##        if '1' in config_type:
##            defaults.append(RELAX_STEP_1)
##        elif '2' in config_type:
##            defaults.append(RELAX_STEP_2)
##    elif 'simulate' in config_type:
##        defaults.append(SIMULATE)
##
##    defaults.append(user_config)
##    return combine(*defaults)

def get_configs(user_configs, config_types):
    """ returns list of complete config dictionaries for simulation / relaxation 
        user_configs is a config dict (or list of config dicts)
        config_types is a string (or list of strings) specifying which default params 
            to apply to the user_configs.
    """
    
    #if isinstance(user_configs, dict):
    #    user_configs = [configs]
    
    if isinstance(config_types, str):
        config_types = [config_types]

    assert isinstance(user_configs, list), "configs must be a list"
    assert isinstance(config_types, list), "config_types must be a list"
    assert len(user_configs) == len(config_types), "user configs and config types not equal"
    assert len(config_types) != 0, "config_types is an empty list"

    return [get_config(*x) for x in zip(user_configs, config_types)]


def get_simulation_configs(user_config):
    """ user_config is a dictionary of oxDNA configuration inputs """
    if isinstance(user_config, list):
        return get_configs(user_config, ['simulate'] * len(user_config))
    else:
        return get_configs(user_config, 'simulate')


def get_relax_configs(step1, step2):
    """ step1 is a dictionary of oxDNA configuration inputs for step 1 of relaxation
        step2 is a dictionary of oxDNA configuration inputs for step 2 of relaxation
    """
    return get_configs([step1, step2], ['relax1', 'relax2'])






class Metrics:
    """ compute metrics based on a .conf file """
    def __init__(self, power, *args):
        """ args is a list of (idx_1, idx_2) tuples of pairs of bases to optimize distance """
        self.base_pairs = args
        self.power = power

    def dist(self, sim_state):
        self.lines = read_conf(sim_state.conf_file)

        return sum([ self.base_dist(*pair)**self.power for pair in self.base_pairs ])

    def base_dist(self, base1_idx, base2_idx):
        """ compute the distance of nucleotides in given conf file
            lines is raw text lines from conf file
        """
        pos1 = get_pos(self.lines[int(base1_idx)])
        pos2 = get_pos(self.lines[int(base2_idx)])

        base_position1 = vector_add(pos1[0], scale(pos1[1], 0.4))
        base_position2 = vector_add(pos2[0], scale(pos2[1], 0.4))
        
        # define POS_BASE 0.4f

        # dr = n1.cm_pos - n2.cm_pos
        # # dr = n1.pos_stack - n2.pos_stack
        # # dr = n1.pos_base - n2.pos_base

        return euclidean_dist(base_position1, base_position2)


def read_conf(conf_path):
    with open(conf_path) as f:
        f.readline()
        f.readline()
        f.readline()
        return f.readlines()


def scale(a, scalar):
    return [x*scalar for x in a]

def euclidean_dist(a, b):
    return sqrt(sum([ x**2 for x in vector_add( a, vector_invert(b) ) ]))

def vector_invert(a):
    return [-x for x in a]

def vector_add(a, b):
    return [a1 + b1 for a1, b1 in zip(a,b)]

def get_pos(line):
    data = [float(x) for x in line.split()]

    # center of mass
    # a1 -- Unit vector indicating orientation of backbone with respect to base
    # a3 -- Unit vector indicating orientation (tilting) of base with respect to backbone
    # velocity
    # angular velocity
    return data[:3], data[3:6], data[6:9], data[9:12], data[12:15]


# if __name__ == '__main__':
#     m = Metrics(sys.argv[1])
#     print(m.dist(*sys.argv[2:]))
