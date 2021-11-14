""" this runs a simulation where tiles are coaxed into binding using some
    sort of optimization method. The program runs, but these methods do not work
"""

import os
import sys
import random
import logging
from pyoxdna import *
from utils import *
from datetime import datetime


logger = None

def get_logger(output_dir):
    global logger

    logger = logging.getLogger('tile_binding')
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(name)s.py %(message)s')
    fh = logging.FileHandler(os.path.join(output_dir, 'search.log'))
    fh.setFormatter(formatter)
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)

    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    sh.setLevel(logging.DEBUG)
    logger.addHandler(sh)

def get_configs(top_file, conf_file, search_steps, gpu=False):
    """ returns relax_step1, relax_step2, test_config, simulation_config """

    step1 = {
        'conf_file': conf_file,
        'topology': top_file,
        'steps': '2000',
    }
    
    step2 = {
        'topology': top_file,
        'print_conf_interval': '50',
    }

    sim_config ={
        ####  GENERIC OPTIONS  ####
        'sim_type': 'MD',
        # 'backend': 'CUDA',
        # 'backend': 'CPU',
        'backend_precision': 'double',
        'restart_step_counter': 0,

        ####  MD OPTIONS  ####
        'dt': '0.005',
        'newtonian_steps': 103,
        'diff_coeff': 2.50,
        #pt = 0.1
        'thermostat': 'john',

        ####  DNA2 OPTIONS ####
        # how salt affects DNA interacitons https://en.wikipedia.org/wiki/DNA#Base_pairing
        # from what I understand, more salt makes it easier for DNA to bind
        # 'interaction_type': 'DNA2',
        # 'salt_concentration': '2',

        ####  SIMULATION OPTIONS  ####
        'steps': search_steps,
        'seed': random.randint(-2147483647, 2147483647), # comment this out for random seed
        'use_average_seq': 0,
        'seq_dep_file': os.path.join(os.environ['OXDNA_HOME'], 'oxDNA1_sequence_dependent_parameters.txt'), 
        'external_forces': 0,
        'T': '45C',
    
        
        # TO TRY
        # 'small_system': 1,
        # 'preserve_topology':1,
        'verlet_skin': .25, # maybe lower this to 1, 0.15 or 0.05
        

        ####  MONTE CARLO OPTIONS  ####
        'ensemble': 'NVT',
        'delta_translation': 0.11, #0.05
        'delta_rotation': 0.05,

        ####  UMBRELLA SAMPLING  ####
        # 'umbrella_sampling': 1,
        # 'op_file': 'op.txt',
        # 'weights_file': 'wfile.txt',


        ####  SKEPTICAL OF THESE  ####
        # 'debug': '1',
        # 'refresh_vel': '1',
        # 'small_system': 1,

        ####    INPUT / OUTPUT    ####
        'conf_file': conf_file,
        'topology': top_file,
        'trajectory_file': 'trajectory.dat',
        'log_file': 'log.dat',
        'energy_file': 'energy.dat',
        'time_scale': 'linear',
        'print_conf_interval': search_steps,
        
        'print_energy_every': search_steps*2, # don't print unnecessary output
        'no_stdout_energy':1,

        'max_io': 20,
    }
    sim_config['backend'] = 'CUDA' if gpu else 'CPU'

    test_config = sim_config.copy()
    test_config['steps'] = 0

    return step1, step2, test_config, sim_config

def log(step, best_conf, dist, seed):
    logger.info(f'step {step}: {os.path.basename(best_conf)} with seed {seed}\t{dist}')

def search(sim, config, name):
    """ runs a random simulation
        sim is pyoxdna object, config is configuration dict, seed are integers
    """
    
    outfile = f'{name}.conf'
    sim.run_one(config, f'{name}.params', outfile)
    return outfile

def get_files(name, output_dir):
    input_file = os.path.join(output_dir, f'{name}.params')
    output_file = os.path.join(output_dir, f'{name}.conf')
    return input_file, output_file

class SimState:
    def __init__(self, sim, config, input_file=None, conf_file=None):
        self.simulation = sim
        self.config = config
        self.input_file = input_file
        self.conf_file = conf_file

    def search(self):

        # get filenames that don't exist yet
        name = 'start'
        input_file, output_file = get_files(name, self.simulation.output_dir)
        while os.path.isfile(input_file) or os.path.isfile(output_file):
            name = str(random.randint(0, 1000000))
            input_file, output_file = get_files(name, self.simulation.output_dir)
            
        config = self.config.copy()
        config['seed'] = get_seed()

        self.simulation.run_one(config, input_file, output_file)

        return SimState(self.simulation, config, input_file, output_file)

    def clean(self):
        for f in [self.conf_file, self.input_file]:
            try:
                rm(f)
            except:
                raise

def clean_dir(directory, current_step, best_conf):
    """ removes all non-best search files from directory """
    current_step = str(current_step)
    best = os.path.splitext(os.path.basename(best_conf))[0]

    for f in os.listdir(directory):
        try:
            name = os.path.splitext(f)[0]
            # remove all files from current step that aren't the best
            if name.split('_')[0] == current_step and name != best:
                rm(os.path.join(directory, f))
        except ValueError:
            continue
        
def relax(output_dir, step1, step2, test_config):
    """ relax a configuration, cleanup when done, and return path to relaxed configuration"""
    relaxation_dir = os.path.join(output_dir, 'relaxation')
    r = DNARelaxer(output_dir=relaxation_dir)
    r.relax(step1, step2, test_config)

    relaxed_conf = os.path.join(output_dir, 'start.conf')
    copy_file(os.path.join(relaxation_dir, 'final.conf'), relaxed_conf)

    rm(relaxation_dir)
    return relaxed_conf
    
def get_seed():
    return random.randint(-2147483647, 2147483647)

def main(output_dir, conf_file, top_file, steps, searches, search_steps, gpu):
    get_logger(output_dir) # gives access to global logger variable
    logger.info(f'using {os.path.basename(conf_file)} for {steps} steps, {searches} searches, {search_steps} search_steps')
    
    # get configs for relaxation and simulation
    step1, step2, test_config, sim_config = get_configs(top_file, conf_file, search_steps, gpu)

    sim_config['conf_file'] = relax(output_dir, step1, step2, test_config)
    
     
    # for cubic_6_tiles
    bases = [
        # glue 1
        [3438,3492],
        [3437,3493],
        [3436,3494],
        [3435,3495],
        [3434,3496],

        # glue 2
        # [3013,3539],
        # [3012,3540],
        # [3011,3541],
        # [3010,3542],
        # [3009,3543],
    ]

    # bases to optimize for small_barrel
    # bases = [
    #     # glue 1
    #     [3177,3595],
    #     # [3178,3594],
    #     # [3179,3593],
    #     # [3180,3592],
    #     # [3181,3591],

    #     # # glue 2
    #     # [3102,3544],
    #     # [3101,3545],
    #     # [3100,3546],
    #     # [3099,3547],
    #     # [3098,3548],
    # ]

    metrics = Metrics(0.5, *bases)
    # metrics = Metrics(1, *bases)
    # metrics = Metrics(2, *bases)

    # create a simulation object
    sim = pyoxdna(output_dir=output_dir)

    # simple_optimize(sim, metrics, sim_config, steps, searches, search_steps)
    minmax_optimize(sim, metrics, sim_config, steps, searches, search_steps)
    
    rm(os.path.join(output_dir, 'energy.dat'))
    rm(os.path.join(output_dir, 'trajectory.dat'))


def minmax_optimize(sim, metrics, sim_config, steps, searches, search_steps):
    """ sim_config['conf_file'] is guaranteed to be a relaxed starting configuration file """

    current_state = SimState(sim, sim_config)
    # children = None
    start_conf = sim_config['conf_file']
    logger.info(f'start: {start_conf}')

    # min_dist = metrics.dist(best_conf)
    # min_dist = None
    # best_seed = None
    
    for i in range(steps):
        
        min_max = None
        best_state = None

        # get the minimum of the maximum child distances
        for j in range(searches):
            
            # filename = os.path.join(sim.output_dir, f'{step}_{search}')
            # conf = search(sim, sim_config, filename)


            new_state = current_state.search()
            
            max_dist = 0
            # get max child distance
            for k in range(searches):
                child = new_state.search()
                dist = metrics.dist(child)
                
                child.clean()
                
                logger.debug(f'{child.conf_file} distance {dist}')

                if dist > max_dist:
                    max_dist = dist
                if min_max is not None and dist > min_max:
                    break
                    
                
            logger.debug(f'{new_state.conf_file} max_dist {max_dist}')

            if min_max is None or max_dist < min_max:
                min_max = max_dist
                best_state = new_state
            else:
                new_state.clean()

        logger.info(f'step {i}: {best_state.input_file} {best_state.conf_file} {min_max}')
        current_state = best_state
        # children = new_children


    # clean_dir(sim.output_dir, step_num, best_conf)
    
    # write step to log file
    # log(step_num, best_conf, min_dist, best_seed)

def simple_optimize(sim, metrics, sim_config, steps, searches, search_steps):
    """ sim_config['conf_file'] is guaranteed to be a relaxed starting configuration file """

    # start with output of relaxation as the best conf so far
    best_conf = sim_config['conf_file']
    # min_dist = metrics.dist(best_conf)
    min_dist = None
    best_seed = None
    
    for step_num in range(steps):
        
        # start step with best .conf found so far
        sim_config['conf_file'] = best_conf
        for search_num in range(searches):
            
            seed = get_seed()

            sim_config['seed'] = seed
            filename = os.path.join(sim.output_dir, f'{step}_{search}')
            conf = search(sim, sim_config, filename)
            
            # calculate optimization distance
            dist = metrics.dist(conf)
            logger.debug(f'{conf} distance {dist}')
            
            # update min_dist and best_conf
            if min_dist is None or dist < min_dist:
                min_dist = dist
                best_conf = conf
                best_seed = seed
            else:
                # technically not needed b/c of clean_dir()
                # but might be good for large number of searches (~1000)
                rm(conf)
            
        clean_dir(sim.output_dir, step_num, best_conf)
        
        # write step to log file
        log(step_num, best_conf, min_dist, best_seed)


def launch_self(conf_file, top_file, number):
    """ submits this script to be run on Razor """

    # walltime='72:00:00'
    walltime='3:00:00'

    launcher = JobLauncher(
        email=EMAIL_ADDRESS,
        queue='gpu72',
        nodes=1,
        ppn='6', 
        walltime=walltime
    )


    for i in range(number):
        timestamp = datetime.now().strftime('%s') + str(i)
        sim_name = f'tile_binding_{timestamp}'
        python_file = os.path.join(CODE_HOME, 'tile_binding.py')
        
        launcher.launch_job(
            sim_name=sim_name,
            command=f'python3 {python_file} run {conf_file} {top_file}', 
            working_dir=os.path.join(SIM_HOME, sim_name)
        )

if __name__ == '__main__':
   
    args = sys.argv[1:]
    
    if args[0] == 'run': # run main program
        # local usage:
        #   python tile_binding.py run conf_file top_file [output_dir]
        conf_file, top_file = args[1:3]

        try:
            output_dir = args[3]
            mkdir(output_dir)
        except:
            output_dir = '.'

        main(output_dir, conf_file, top_file, steps=30, searches=10, search_steps=500, gpu=False)
        # main(output_dir, conf_file, top_file, steps=20, searches=5, search_steps=1000, gpu=False)
    
    else: # launch program in pbs script

        # usage: 
        #   python launch_tile_bind.py conf_file top_file number_simulations
        #   cd; ./python3.sh launch_tile_bind.py "/home/krs028/oxDNA-simulations/kyle/my_code/oxdna_files/small_barrel.conf /home/krs028/oxDNA-simulations/kyle/my_code/oxdna_files/small_barrel.top 1"
        conf_file, top_file = [os.path.abspath(x) for x in args[:2]]
        number = int(args[2])

        launch_self(conf_file, top_file, number)


