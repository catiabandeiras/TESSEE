# main.py - the main file where the simulations will be launched

import os

import shutil

import simpy

from stemfactory import *

from databases import *

number_runs = 1

print('Started simulation')

# Starts by deleting the auxiliary directories
# From previous simulation

cwd = os.getcwd()

csv_dir = cwd + '/csv_files'

plot_dir = cwd + '/plots'

shutil.rmtree(csv_dir, ignore_errors = True)
shutil.rmtree(plot_dir, ignore_errors = True)

# Reads the database and starts the factory

for i in range(number_runs):

    env = simpy.Environment()

    db = Database()

    stemfactory = StemFactory(env, db)

    env.process(stemfactory.run(env, db))

    env.run()

print('Ended simulation')
