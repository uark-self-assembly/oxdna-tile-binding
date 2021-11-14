import os
import sys

""" this program converts the output of optimization_binding.py
into a trajectory file which can be viewed in oxdna-viewer """

def combine_confs(working_dir):
    conf_files = ['start.conf'] + sorted([f for f in os.listdir(working_dir) if '.conf' in f and '_' in f], key=lambda x: int(x.split('_')[0]))
    print(conf_files)
    with open('compiled_trajectory.dat', 'w+') as f:
        for conf in conf_files:
            with open(os.path.join(working_dir, conf), 'r') as c:
                f.write(c.read())


if __name__ == '__main__':
    assert len(sys.argv) == 2, f"""usage: python {sys.argv[0]} conf_dir
conf_dir is a directory with configuration files to be combined
the starting conf file is named start.conf and the rest are named <index>_<name>.conf
non-starting conf files are ordered by index.
"""
    working_dir = sys.argv[1]
    combine_confs(working_dir)
