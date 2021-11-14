"""
    relaxes and runs a simulation to quantify the computation profile of oxDNA on dford's GPUs
    the average computation time of num_simulations simulations is written
        in <working_dir>/results.txt and ~/results.csv
    usage: python run.py output_dir num_strands
    
"""


import os
import sys
import random
from utils import *
from time import perf_counter, process_time
from datetime import datetime
from pyoxdna import DNARelaxer, pyoxdna
from pyoxdna.utils import current_time, copy_file, mkdir



def create_files(output_dir, num_strands, bases_per_strand, box_size):
    """ creates .top and .conf files. returns the filepath prefix for the files 
        basename of prefix is the number of base pairs
    """
    top_file = os.path.join(output_dir, 'sim.top')
    conf_file = os.path.join(output_dir, 'sim.conf')

    with open(top_file, 'w') as f:
        f.write(get_top_file(num_strands, bases_per_strand))

    with open(conf_file, 'w') as f:
        f.write(get_conf_file(num_strands, bases_per_strand, box_size))

    return top_file, conf_file

def get_top_file(num_strands, bases_per_strand):
    """ returns the contents of the top file with num_strands """

    output = str(num_strands*bases_per_strand) + ' ' + str(num_strands) + '\n'

    for i in range(num_strands):
        for j in range(bases_per_strand):
            nxt = i*bases_per_strand+j+1 if j < bases_per_strand - 1 else -1
            prev = i*bases_per_strand+j-1 if j != 0 else -1
            output += str(i+1) + ' ' + random_base() + ' ' +str(prev)+ ' '+str(nxt) + '\n'

    return output

def get_conf_file(num_strands, bases_per_strand, box_size):
    """ returns the contents of the conf file with num_strands """
    
    output = 't = 0\nb = ' + str(box_size) + '.00 ' + str(box_size) + '.00 ' + str(box_size) + '.00\nE = 0.00 0.00 0.00\n'

    position = generate_xyz()
    joined_spacing = 0.526329901992871 #for joined nucleotides
    z_spacing = 2
    # close_spacing = 0.60 #for joined nucleotides
    spacing = bases_per_strand*joined_spacing + 2 # for vertically stacked strands
    # spacing = 1

    for i in range(num_strands):
        center = [value for value in next(position)]
        for j in range(bases_per_strand):
            # output += ' '.join(['{:.2f}'.format(x) for x in [center[0]+j*close_spacing, center[1], center[2]]]) \
            output += ' '.join(['{:.2f}'.format(x + box_size/2) for x in [center[0]*spacing+j*joined_spacing, center[1]*z_spacing, center[2]]]) \
                + ' {:.2f} {:.2f} {:.2f}'.format(0, 1, 0) \
                + ' {:.2f} {:.2f} {:.2f}'.format(0, 0, 1) \
                + ' {:.2f} {:.2f} {:.2f}'.format(0, 0, 0) \
                + ' {:.2f} {:.2f} {:.2f}'.format(0, 0, 0)+'\n'

    return output

def generate_xyz():
    # make a 3D "spiral" :)
    # this is technically not a spiral, but a 3D cube growing pattern
    discrete_dist = 0

    while True:
        points = generate_points(discrete_dist)
        point = next(points, None)

        while point is not None:
            yield point
            point = next(points, None)

        discrete_dist += 1

def generate_points(discrete_dist):

    # faces of the cube (one coordinate is at an extreme)
    for i in range(1-discrete_dist, discrete_dist):
        for j in range(1-discrete_dist, discrete_dist):
            yield [i, j, discrete_dist]
            yield [i, j, -discrete_dist]
            
            yield [j, discrete_dist, i]
            yield [j, -discrete_dist, i]

            yield [-discrete_dist, i, j]
            yield [discrete_dist, i, j]
    
    # edges of the cube (two coordinates are at extremes)
    for i in range(1-discrete_dist, discrete_dist):
        yield [discrete_dist, discrete_dist, i]
        yield [-discrete_dist, discrete_dist, i]
        yield [discrete_dist, -discrete_dist, i]
        yield [-discrete_dist, -discrete_dist, i]

        yield [discrete_dist, i, discrete_dist]
        yield [-discrete_dist, i, discrete_dist]
        yield [discrete_dist, i, -discrete_dist]
        yield [-discrete_dist, i, -discrete_dist]

        yield [i, discrete_dist, discrete_dist]
        yield [i, -discrete_dist, discrete_dist]
        yield [i, discrete_dist, -discrete_dist]
        yield [i, -discrete_dist, -discrete_dist]

    # vetices of the cube (all coordinates are at extremes)
    if discrete_dist == 0:
        yield [0, 0, 0]

    else:
        yield [discrete_dist, discrete_dist, discrete_dist]
        yield [-discrete_dist, -discrete_dist, -discrete_dist]

        yield [discrete_dist, discrete_dist, -discrete_dist]
        yield [-discrete_dist, -discrete_dist, discrete_dist]

        yield [discrete_dist, -discrete_dist, discrete_dist]
        yield [-discrete_dist, discrete_dist, -discrete_dist]

        yield [discrete_dist, -discrete_dist, -discrete_dist]
        yield [-discrete_dist, discrete_dist, discrete_dist]

def random_base():
    return random.choice(['A', 'T', 'G', 'C'])



class ResultsWriter:
    def __init__(self, csv_file):
        self.csv_file = csv_file

    def write_results(self, num_strands, sim_seconds_total, sim_seconds_process, gpu):
        with open(self.csv_file, 'a+') as f:
            f.write(str(num_strands)+','+str(sim_seconds_total)+','+str(sim_seconds_process)+','+str(gpu)+','+str(datetime.now())+'\n')

        print(num_strands, sim_seconds_total, sim_seconds_process, gpu, datetime.now())

        with open('results.txt', 'w') as f:
            f.write(str(num_strands)+','+str(sim_seconds_total)+','+str(sim_seconds_process)+','+str(gpu)+','+str(datetime.now())+'\n')



def get_configs(top_file, conf_file, gpu=False):
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
        'steps': 5e3,
        # 'seed': -694528013, # comment this out for random seed
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
        'print_conf_interval': 1e10,
        'print_energy_every': 1e10,
        'no_stdout_energy':1,

        'max_io': 20,
    }

    test_config = sim_config.copy()
    test_config['steps'] = 0

    # we want to relax with the cpu
    test_config['backend'] = 'CPU'

    sim_config['backend'] = 'CUDA' if gpu else 'CPU'

    return step1, step2, test_config, sim_config

def time_simulations(sim, sim_config, num_simulations):
    time_sum_total = 0
    time_sum_process = 0
    for i in range(num_simulations):
        start_total = perf_counter()
        start_process = process_time()
        sim.run(sim_config)
        time_sum_total += perf_counter() - start_total
        time_sum_process += process_time() - start_process
    return time_sum_total/num_simulations, time_sum_process

def main(output_dir, num_strands, bases_per_strand, num_simulations, box_size, gpu):
    assert (box_size / 5)**3 >= num_strands
    
    mkdir(output_dir)
    relaxation_dir = os.path.join(output_dir, 'relaxation')

    # create sim.conf and sim.top files
    top_file, conf_file = create_files(output_dir, num_strands, bases_per_strand, box_size)
    
    # get configs for relaxation and simulation
    step1, step2, test_config, sim_config = get_configs(top_file, conf_file, gpu)

    if(not os.path.isdir(relaxation_dir)):
        # relax sim.conf
        print("Relaxing...\n")
        r = DNARelaxer(output_dir=relaxation_dir)
        r.relax(step1, step2, test_config)
        print("Done relaxing. Begining simulation\n")
    
    # add output of relaxation to sim configuration file
    sim_config['conf_file'] = os.path.join(relaxation_dir, 'final.conf')

    # create a simulation and measure the average time per simulation 
    # sim = pyoxdna(output_dir=output_dir, output_level=0)
    sim = pyoxdna(output_dir=output_dir, output_level=1)
    avg_time_total, avg_time_process = time_simulations(sim, sim_config, num_simulations)
    print("Simulation complete, writing results")

    # write results to ~/results.csv, print them, and write them in results.txt
    rw = ResultsWriter(os.path.join(os.environ['HOME'], 'results.csv'))
    rw.write_results(num_strands, avg_time_total, avg_time_process, gpu)


def launch_self(start, end):
    walltime='72:00:00'

    launcher = JobLauncher(
        email=EMAIL_ADDRESS,
        queue='gpu72',
        nodes=1,
        ppn='6', 
        walltime=walltime
    )

    for i in range(start, end):
        num_strands = 2**i
        
        sim_name = f'simulation_{num_strands}'
        current_file = os.path.abspath(__file__)

        working_dir=os.path.join(SIM_HOME, sim_name)
        
        launcher.launch_job(
            sim_name=sim_name,
            working_dir=working_dir,
            command=f'python3 {current_file} run {num_strands} {working_dir}'
        )
    


if __name__ == '__main__':
    # seg fault at 1 x 1, 4 x 734 and 10 x 32 bases
    # doesn't seg fault with bases_per_strand = 2

    args = sys.argv[1:]

    if args[0] == 'run':
        # usage python experiment.py run num_strands [output_dir]
        box_size = 1000
        gpu = True
        bases_per_strand = 2
        num_strands = int(args[1])
        
        try:
            output_dir = args[2]
        except:
            output_dir = '.'

        if num_strands == 16384:
            num_simulations = 200
        elif num_strands == 32768:
            num_simulations = 150
        elif num_strands == 65536:
            num_simulations = 100
        elif num_strands == 131072:
            num_simulations = 100
        else:
            num_simulations = 500
        
        main(output_dir, num_strands, bases_per_strand, num_simulations, box_size, gpu)

    else:
        # usage: python3 launch_experiment.py start end
        launch_self(*[int(x) for x in args[:2]])
