'''main.py - the main file where thelif self.passage_no == self.max_passages_exp:e simulations will be launched'''

import simpy

from stemfactory import *

#from databases_esc import *

from databases_ipsc import *

number_runs = 1

print('Pluripotent stem cells, Semma protocol, obtaining islet cells, planar, regular, simulation started, allogeneic')

for i in range(number_runs):

	env = simpy.Environment()

	#Create the initial database

	db = Database()

	stemfactory = StemFactory(env,db)

	env.process(stemfactory.run(env,db))

	env.run()

print('Pluripotent stem cells, Semma protocol, obtaining islet cells, planar, regular, allogeneic')