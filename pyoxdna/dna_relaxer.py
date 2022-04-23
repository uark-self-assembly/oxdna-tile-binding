import os
import sys
import signal
from .utils import *
from .pyoxdna import pyoxdna
from traceback import print_exc

def get_file(path, num):
    filename, ext = os.path.splitext(path)
    return filename + '-' + str(num) + ext

class DNARelaxer:
    """ class for relaxing a DNA simulation until it passes a test simulation """

    def __init__(self, output_dir='output', exit_on_passing=True, autoclean=True):
        """ @input output_dir directory where input / output files are stored, defaults to 'output'
            @input name name of simulation, defaults to 'relax'
            @input exit_on_passing if True, the program exits once a relaxed.conf passes self._test
                (i.e. the output is the first passing relaxed configuation found)
        """
        self.autoclean = autoclean
        self.output_dir = output_dir
        self.child_processes = []
        self.exit_on_passing = exit_on_passing
        self.pids_to_kill = os.path.join(self.output_dir, '.pids_to_kill'+current_time())

    def relax(self, rstep1, rstep2, testing_config, max_attempts=-1):
        """ @input rstep1 config dictionary for relaxation step 1
            @input rstep2 config dictionary for relaxation step 2
            @input testing_config config dictionary for testing relaxation outputs
            @input max_attempts integer to bound number of attempts, default -1 for no bound
        """
        print_message(5, "About to launch relax_process...")
        relax_process = self._launch_subprocess(self._relax, (rstep1, rstep2))
        print_message(5, "... relax_process launched")

        # keep track of when output file changes, copying to new file when it does.
        # test each new output to determine when it is relaxed enough
        last_change = 0
        counter = 0

        output_path = os.path.join(self.output_dir, rstep2['lastconf_file'])

        msg = 'In DNARelaxer, about to begin relaxation loop'
        if max_attempts > -1:
            msg = msg + f' for a maximum of {max_attempts} iterations'
        print_message(3, msg)
        while relax_process.is_alive() and not counter == max_attempts:
            sleep(1)
            
            try:
                timestamp = os.path.getmtime(output_path)
                if last_change != timestamp:
                    new_file = get_file(output_path, counter)
                    copy_file(output_path, new_file)

                    self._run_test(new_file, testing_config)

                    last_change = timestamp
                    counter += 1
                    
            except FileNotFoundError:
                pass # output file has not been created yet
        
        if counter == max_attempts:
            print(f"Relaxation reached the maximum allowed number of attempts, {max_attempts}")
    
        self.exit()
        
    def exit(self):
        for proc in self.child_processes:
            try:
                kill(proc.pid)
            except AttributeError:
                pass # process has already exited
        
        try:
            with open(self.pids_to_kill, 'r') as f:
                kill_pid_file(f.read())
        except FileNotFoundError:
            pass

        rm(self.pids_to_kill)

        
    def _run_test(self, conf_file, testing_config):
        """ tests if the DNA structure is relaxed enough
            if the test passes and self.exit_on_passing is true, kills program
            conf_file is name of input .conf file
            testing_config is dict of oxDNA params
        """
        conf_file = os.path.abspath(conf_file)
        print(f"Testing relaxation of configuration file '{conf_file}'")
        testing_config['conf_file'] = conf_file
        self._launch_subprocess(self._test, (testing_config,))
        
    def _launch_subprocess(self, function, fn_args):
        """ launches a subprocess """
        process = Process(target=function, args=fn_args)
        self.child_processes.append(process)
        process.start() # don't block
        return process

    def _test(self, testing_config):
        """ runs an oxdna simulation according to testing_config oxdna settings
            if the simulation does not error, it is assumedthat the input .conf file is relaxed enough
            quits the program if the simulation doesn't error
        """
        conf_file = testing_config['conf_file']
        output_dir = os.path.join(self.output_dir, 'test_'+current_time())
        sim = pyoxdna(output_dir=output_dir, output_level=1 if self.autoclean else 3)

        try:
            sim.run(testing_config)
            print_message(5, conf_file + ' test simulation PASSED')
            copy_file(conf_file, os.path.join(os.path.dirname(conf_file), 'final.conf'))

            if self.autoclean:
                rm(output_dir)
                rm(conf_file)

            if self.exit_on_passing:
                self.exit()

        except CalledProcessError:
            print_message(5, conf_file + ' test simulation FAILED')
            
            if self.autoclean:
                rm(output_dir)
                rm(conf_file)


    def _relax(self, relax_step_1, relax_step_2):
        """ @input rstep1 config dictionary for relaxation step 1
            @input rstep2 config dictionary for relaxation step 2
        """
        
        sim = pyoxdna(output_dir=self.output_dir)
        with open(self.pids_to_kill, 'w') as f:
            f.write(sim.pid_file)
        sim.relax(relax_step_1, relax_step_2)

    def __del__(self):
        self.exit()
        print('exiting...')
