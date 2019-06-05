# expansiontechnology.py - attributes associated to each individual flask
# Inherits attributes from the donor class

import sys
import simpy
import math
import numpy as np
import pandas as pd
from numpy.random import *
from scipy.integrate import odeint

from databases import *

class ExpansionTechnology(object):

    def __init__(self, name, sd, donor, env, db, number_wells):

        self.name = name  # Identifies the name of the technology

        # Starts by creating the areas of the specific
        # Technologies that can be involved in iPSC culture

        if name == '6well':

            self.area = 19.6*6

        elif name == '500spinner':

            self.area = 500  # Just a dummy variable for now

        else:

            self.area = db.area_planar_ets[name]  # Comes from the database

        self.seeding_density = sd

        # INITIALIZE PARAMETERS

        self.Xv0 = 0

        self.inc_time = donor.inc_time_passage_exp[donor.passage_no]

        self.time_to_confluence = donor.time_to_confluence[donor.passage_no]

        self.final_cells = 0

        self.expanded = 0

        # INTRA-DONOR VARIABILITY
        # hd - harvesting density
        # gr - growth rate
        # If set to 1, parameters do not vary from the donor specific baseline

        self.intra_donor_hd_factor = 1
        self.intra_donor_gr_factor = 1

        self.growth_rate = donor.growth_rate_exp[donor.passage_no]*self.intra_donor_gr_factor

        self.harvest_yield = 0.9  # Reduction of the number of cells due to trypsinization

        # INITIALIZES DIFFERENTIATION STAGE PARAMETERS
        # Only applicable if the process is related to pluripotent
        # stem cells for now

        self.current_diff_stage = 0

        self.number_wells = number_wells

    def expansion(self, donor, env, db):

        # Expand the given ET type

        # Start by providing the initial cell number of this ET
        # self.Xv0

        self.Xv0 = math.floor(self.seeding_density*db.area_planar_ets[self.name])

        # Seeding of the flask, get the volume media spent

        donor.volume_media_spent_expansion += db.media_volume_planar_ets[self.name]

        # Advances the time of simulation
        # By the time it takes to perform a seeding operation 

        yield env.timeout(db.et_operation_times_planar.loc[self.name,'seeding'])

        # In this setting, the flask is active until a specified time for confluence
        # is reached. However, a specific cell number can also be used

        while self.inc_time < donor.time_to_confluence[donor.passage_no]:

            # Advances the time to the next feeding
            # It simulates how the flasks are checked for cell numbers
            # and reagent/metabolite concentration

            yield env.timeout(db.time_bw_feedings)

            self.inc_time += db.time_bw_feedings

            # Accounts for specific growth rates after lag phase

            if self.inc_time > donor.lag_phase_exp:

                # Number of cells will grow according to time
                # and specific growth rates

                self.final_cells = self.Xv0*math.exp(self.growth_rate*(self.inc_time-donor.lag_phase_exp))

            else:

                # If this is still on lag phase
                # Number of cells is kept constant to the initial

                self.final_cells = self.Xv0


            if self.inc_time >= donor.time_to_confluence[donor.passage_no]:

                # Time has passed, so flask is going for harvesting

                break


            else:

                # If needs to advance time
                # Account for consumption of culture media after exchange

                yield env.timeout(db.et_operation_times_planar.loc[self.name,'feeding'])

                donor.volume_media_spent_expansion += db.media_volume_planar_ets[self.name]

                self.inc_time += db.et_operation_times_planar.loc[self.name,'feeding']

        # Harvesting with the reagent (trpysin like) and baseline culture media

        donor.volume_tryple_spent_expansion += db.harvesting_volume_planar_ets[self.name] # mudar para trypLE ou harvesting reagent

        # TrypLE requires a certain amount of media to inactivate it
        # Thus far, assume it is 2 times the volume of harvesting medium

        donor.volume_media_spent_expansion += 2*db.harvesting_volume_planar_ets[self.name]

        yield env.timeout(db.harvesting_time_incubator)  # Advances time to flask processing

        # Collects a certain number of cells after the harvest yield

        cells_after_harvesting = self.final_cells*self.harvest_yield

        # Adds the collected cells to the cells from that donor in the specific passage

        donor.cells_per_passage_exp[donor.passage_no] += cells_after_harvesting

        self.expanded = 1

        donor.harvested_per_passage_exp[donor.passage_no] += 1  # Increases no. flasks per passage

    def growthcurve(self):

        # Method that could allow different growth curves
        # Including Ordinary Differential Equations (CDE)

        miu_max = self.growth_rate*self.intra_donor_gr_factor

        #use the time t as final for curve integration!

        t = np.linspace(0,self.inc_time,num=100)

        cellgrowth_int_handle = lambda Xv,t: self.growthcurve_ode(Xv,t,miu_max)

        cell_per_time = odeint(cellgrowth_int_handle,self.Xv0,t)

        #Returned cell numbers have to be integers

        return math.floor(cell_per_time[-1,-1])

    def growthcurve_ode(self,Xv,t,k):

        # The most simple method for cell growth

        return k*Xv

    def differentiation(self,donor,env,db):

        # Determine what is the initial number of PSCs seeded in each 6 well plate

        if donor.diff_scheme == 'planar':

            # Related to multi well plates

            volume_required = db.media_volume_well * self.number_wells

            volume_tryple = db.tryple_volume_well * self.number_wells

        elif donor.diff_scheme == 'spinner':

            # Assumes a volume based scheme for the spinner

            volume_required = db.media_volume_spinner

            volume_tryple = db.tryple_volume_spinner

        # Inputs the initial cell numbers

        self.Xv0 = math.floor(self.seeding_density*volume_required)

        # Start the feeding and then advance it all.
        # Assume that the volume is always the same for every day that the medium is changed.

        # Specific aggregation media is necessary

        donor.volume_media_spent_aggregation += volume_required*db.time_aggregation

        yield env.timeout(db.time_aggregation)  # Advances time needed for aggregation

        # DIFFERENTIATION STAGE INITIALIZATION

        diff_time = 0

        self.current_diff_stage = 1

        final_diff_time = db.time_stage.values[-1]  # The last day of the differentiation checkpoints

        # Iterate across all stages

        for stage in range(1,db.no_diff_stages + 1):

            # Spends the media from the specific stage
            # Each specific stage media has a specific cost

            donor.volume_media_spent_differentiation[stage-1] += volume_required*db.number_feedings_stage[stage-1]

            # Advances the specific time

            yield env.timeout(db.time_stage[stage-1] - diff_time)

            diff_time += (db.time_stage[stage-1]-diff_time)

            self.current_diff_stage += 1  # Advances the differentiation stages

        # Spends tryple to detach the cells from the coating

        donor.volume_tryple_spent_differentiation += volume_tryple

        donor.harvested_diff += 1

        # Assuming that there is no loss or growth
        # The differented cells are passed with the same number
        # from the aggregation

        donor.cells_per_passage_dif[donor.max_passages_dif-1] += self.Xv0

        self.differentiated = 1


    def diff_stage_selector(self,donor,env,db,diff_time):

        # Select what is the number of the current stage of differentiation

        number_stage = 0

        while number_stage < db.no_diff_stages:

            if diff_time < db.time_stage.values[number_stage]:

                self.current_diff_stage = number_stage + 1

                

                break

            else:

                number_stage += 1

