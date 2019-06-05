# donor.py - class responsible for the coordination of the donors

import sys
import simpy
import random
import math
import logging
import time
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from numpy.random import *
from scipy import integrate, interpolate
from scipy import optimize

# Import the aditional modules

from databases import *

from expansiontechnology import *

class Donor(object):

    def __init__(self, env, donor_index, db):

        self.env = env # Necessary parameter for discrete event simulation

        self.donor_index = donor_index

        # INITIALIZING UNIT OPERATIONS PARAMETERS

        self.passage_no = 0

        self.isolated = 0

        self.expanded = 0

        self.finished = 0

        self.max_passages_exp = db.max_passages_exp

        self.max_passages_dif = 1 
        
        self.doses_per_donor = db.doses_per_batch  # Maximum number of doses

        self.cells_per_dose = db.cells_per_dose

        # Record the different cells of the passage

        self.cells_per_passage_exp = [0]*(self.max_passages_exp+1)
        self.cells_per_passage_dif = [0]*(self.max_passages_dif+1)

        # SEEDING DENSITIES

        self.base_sd_exp = 3000  # cells/cm^2
        self.base_sd_diff = 1e6  # cells/ml

        # Initialize the simulations with the growth rate and base densities at the same as expansion

        self.base_sd = self.base_sd_exp

        # GROWTH DATES (per day)

        self.growth_rate_exp = [0.3]*(self.max_passages_exp+1)  # Example values
        self.growth_rate_dif = [0.3]*(self.max_passages_exp+1)  # Example values

        # TIMES TO CONFLUENCE

        # For now, assume it's only based on the expansion protocols

        self.lag_phase_exp = 3

        self.time_to_confluence = db.confluence_passage_exp

        self.Xv0 = 0  # Initial number of cells

        # INITIALIZE AUXILIARY STORAGE ARRAYS

        self.inc_time_passage_exp = [0]*(self.max_passages_exp+1)

        self.harvested_per_passage_exp = [0]*(self.max_passages_exp+1)

        self.et_list = []

        self.processed = 0

        self.cryovials_stored = 0

        self.tested = 0

        self.last_passage_name = ''

        self.last_passage_number = 0

        # INITIALIZE CONSUMABLES EXPENDITURES

        self.planar_ets_spent = pd.Series(data = [0]*len(db.names_of_planar_ets),
            index = db.names_of_planar_ets)

        self.planar_ets_spent_isolation = pd.Series(data = [0]*len(db.names_of_planar_ets),
            index = db.names_of_planar_ets)


        self.planar_ets_spent_expansion = pd.Series(data = [0]*len(db.names_of_planar_ets),
            index = db.names_of_planar_ets)

        self.cryovials_spent = 0


        # INITIALIZE REAGENT EXPENDITURES

        self.volume_media_spent = 0

        self.volume_hr_spent = 0 # hr = harvesting reagent

        self.volume_media_spent_isolation = 0

        self.volume_tryple_spent_isolation = 0

        self.volume_media_spent_expansion = 0

        self.volume_tryple_spent_expansion = 0

        self.volume_media_spent_DSP = 0

        self.volume_FP_spent = 0  # FP = Ficoll-Paque for density gradient centrifugation

        self.volume_PBS_spent = 0  # PBS = buffer

        self.cryomedium_spent = 0

        # INITIALIZE TIMES

        self.initial_time_processing = 0

        self.end_isolation_time = 0

        self.end_expansion_time = 0

        self.end_differentiation_time = 0

        self.end_DSP_time = 0

        self.time_isolation = 0

        self.time_expansion = 0

        self.time_differentiation = 0

        self.time_processing = 0

        self.time_DSP = 0

        # INITIALIZE COSTS

        self.total_costs_reprogramming = 0 # It's 0 for ESCs expansion, but some value when buying a vial of iPSCs

        self.total_costs_isolation = 0

        self.total_costs_expansion = 0

        self.total_costs_differentiation = 0

        self.total_costs_DSP = 0

        self.total_costs_release = 0

        self.total_costs_ets = 0

        self.total_costs_reagents = 0

        self.total_costs_donor = 0

        self.total_hpl_costs = 0

        # INITIALIZE CELLS

        initial_pscs = 1e6 #Make a dummy variable for either initial ESCs or the iPSCs coming from a vial

        self.cells_P0 = initial_pscs

        self.cells_per_passage_exp[0] = self.cells_P0

        self.doses_no = 0

    def processing(self,env,db):

        # Method that coordinates the complete processing of each donor

        self.initial_time_processing = env.now

        #print('Initial processing time %d' % self.initial_time_processing)

        # Assuming that reprogramming is a service, in case it happens, launch the expansion process

        exp = self.expansion(env,db)

        env.process(exp)

        while self.expanded == 0:

            yield env.timeout(0.0001)

        ##Include the differentiation protocol
        ##If there is no need for differentiation (like this adult case), make a turnaround

        # dif = self.differentiation(env,db) #add the process

        # env.process(dif)

        # while self.differentiated == 0: #add the differentiation potential

        #   yield env.timeout(0.0001)

        # # # # yield env.timeout(self.exp_time)

        self.differentiated = 1

        yield env.timeout(0.0001)

        # Launch the downstream process

        dsp = self.downstream(env,db)

        env.process(dsp)

        while self.finished == 0:

            yield env.timeout(0.0001)


        # Launch the quality controls

        qcr = self.release(env,db)

        env.process(qcr)

        while self.tested == 0:

            yield env.timeout(0.0001)

        # Concluded donor processing and adds doses


        print('Processing of donor %d is finished after %d days' % (self.donor_index,env.now))

        self.time_processing = self.time_expansion + self.time_differentiation + self.time_DSP + db.qc_time

        # Calculate the costs

        self.cost_calculator_phase(env, db)

        self.processed = 1

        db.processed_donors += 1

        db.total_doses += self.doses_per_donor


    def expansion(self,env,db):

        #Initialize the passage number and number of doses produced

        self.passage_no = 1

        doses_no = 0

        # Selects for new consumables until the number of passages and doses are reached

        while self.passage_no <= self.max_passages_exp and doses_no <= self.doses_per_donor:

            # Update the initial cells of the actual passage

            # Create the consumables with expansion requirements and seeding densities

            self.create_ets(env,db,'expansion')

            for et in self.et_list:

                # Launches the expansion process of each flask

                et_exp = et.expansion(self,env,db)

                env.process(et_exp)

                yield env.timeout(0)

            while self.harvested_per_passage_exp[self.passage_no] < len(self.et_list):

                # Holds the process off until all the simultaneous flasks per passage are spent

                yield env.timeout(0.0001)


            # Calculate the associated CPDs

            self.cpds = round(math.log(self.cells_per_dose*self.doses_per_donor/self.initial_cells_for_expansion)/math.log(2),2)

            # Stops the simulation if the final passage number is reached

            if self.passage_no == self.max_passages_exp:

                break

            else:

                self.passage_no += 1

        self.expanded = 1  # Finished the expansion

        self.end_expansion_time = env.now

        self.time_expansion = self.end_expansion_time - self.end_isolation_time

    def differentiation(self,env,db):

        # Method that can be adapted for differentiation of stem cells into several types

        self.create_ets(env,db,'differentiation')

        initial_pscs = self.cells_per_passage[self.passage_no]  # Starts with the cells from expansion

        differentiation_time = 0

        differentiation_time_end = db.time_stage[end]

        # Typical differentiation scheme using well plates

        number_wells = initial_pscs/(db.media_volume_well * self.base_sd_diff)

        number_plates = number_wells / 6

        # If this is a process for adult cells, the timeout is automatic

        yield env.timeout(0.001)

    def downstream(self,env,db):

        #For simplicity, use the centrifuge as common just as a proof of concept
        #Use the fixed VR time

        washes_done = 0

        # Perform several washes

        while washes_done < db.number_washes:

            self.volume_media_spent_DSP += self.last_passage_number*db.media_volume_planar_ets[self.last_passage_name] #Needs to be corrected!

            yield env.timeout(db.wash_vr_time)

            washes_done += 1

        # Calculate the number of cryovials for storage 

        cryovials_no = math.ceil(self.cells_per_passage_exp[self.passage_no]/(db.cell_concentration*db.cryovial_vol))

        # Updates process consumables and reagents expenditures

        self.cryovials_spent += cryovials_no

        self.volume_media_spent_DSP += cryovials_no*db.cryovial_vol*(1-db.cryomedium_ratio)

        self.cryomedium_spent += cryovials_no*db.cryovial_vol*db.cryomedium_ratio

        self.cryovials_stored = cryovials_no

        yield env.timeout(db.finish_time)

        #print('Cells processed from donor %d after DSP are %d' % (self.donor_index,self.cells_per_passage[self.passage_no]))

        self.end_DSP_time = env.now

        self.time_DSP = self.end_DSP_time - self.end_expansion_time

        self.finished = 1


    def release(self,env,db):

        # Process related to passing the release or not

        self.cryovials_stored -= 1  # Sacrifice one of the vials for testing

        yield env.timeout(db.qc_time)

        # Draw from a binomial distribution if it's a pass or a fail

        switch = np.random.binomial(1,db.pass_release_ratio)

        if switch == 1:

            print('Release testing of donor %d passed.' % self.donor_index)

            self.doses_no = math.floor(self.cryovials_stored * db.cell_concentration * db.cryovial_vol /self.cells_per_dose)

            print('Doses of donor %d are %d' % (self.donor_index,self.doses_no))

            self.doses_per_donor = self.doses_no

            print('Processing of donor %d finished at passage %d with %d doses produced' % (self.donor_index,self.passage_no,self.doses_no))

            print('Cells from donor %d are stored in %d cryovials with a concentration of %.2f cells/ml' % (self.donor_index,self.cryovials_stored,db.cell_concentration))

        else:

            print('Release testing of donor %d failed.' % self.donor_index)

            # Production of doses is cancelled. All goes to disposal!

            self.doses_per_donor = 0

        self.tested = 1

    def create_ets(self,env,db,stage):

        # First of all, choose what is the parametrization

        if stage == 'expansion':

            # Collects the cells to be divided by the most approximate area of flasks
            # Aims at minimizing incubator occupation to save resources

            self.base_sd = self.base_sd_exp

            cells_to_seed = self.cells_per_passage_exp[self.passage_no-1]

            cells_per_et = self.base_sd*db.area_planar_ets

            no_ets = np.floor(cells_to_seed/cells_per_et)

            #Replace the infinite values by zero

            #print('Number of expansion technologies generated')
            #print(no_ets)

            #Choose the index

            et_name = ''
            sd_error = 1e16
            sd_new = 0

            #Choose, from options that yield the minimum expansion technologies, the ones that waste less cells

            #Subset first so that only values larger than zero appear. Sometimes the best option is not to have a zero

            no_ets_subset = no_ets[no_ets>0]

            minimum_ets = min(no_ets_subset) #So that the minimum number of ETs is never zero!

            for name in db.names_of_planar_ets:

                if no_ets[name] == minimum_ets:

                    sd_per_et = cells_to_seed/(db.area_planar_ets[name]*minimum_ets)

                    sd_error_aux = abs(sd_per_et-self.base_sd)

                    if sd_error_aux < sd_error:

                        sd_error = sd_error_aux

                        et_name = name

                        sd_new = self.base_sd  # Since no relationships between seeding density are incorporated, we will keep the seeding density


            # Creates a limitation in the number of flasks if the total incubator spaces are full

            if stage == 'expansion' and no_ets[et_name] > db.units_resource_planar.loc[et_name,'incubator']:

                no_ets[et_name] = db.units_resource_planar.loc[et_name,'incubator']

            # Creates a list of flasks of the same type by calling the ExpansionTechnology class
            # The components of the list queue up for room in the facility

            # Since this is expansion, the number_wells arg is set to 0

            self.et_list = [ExpansionTechnology(et_name, sd_new, self, env, db, 0) for index in range(int(no_ets[et_name]))]

            # For each passage, the number of flasks generated are added into the database

            self.planar_ets_spent_expansion[et_name] += int(no_ets[et_name])

            # Assume that this is the last passage until further information
            # The last passage number if updated dynamically
            # Possibility that the passages stop earlier for reaching the total demand

            self.last_passage_name = et_name

            self.last_passage_number = no_ets[et_name]

            # At the initial passage, there might be some cell waste/storage to comply
            # with keeping the seeding density constant
            # Models on cell growth rates vs seeding density can be developed
            # given that enough data are made available

            if self.passage_no == 1:

                self.initial_cells_for_expansion = round(no_ets[et_name]*db.area_planar_ets[et_name]*sd_new,0)

                print('Seeding density at P1 %d cells/cm^2' %sd_new)
                print('Initial cells for expansion seeded at P1 are %d' % self.initial_cells_for_expansion)


        elif stage == 'differentiation':

            # For differentiation, only one "passage" is required

            self.base_sd = self.base_sd_diff

            cells_to_seed = self.cells_per_passage_exp[self.passage_no]

            

    def cost_calculator_phase(self,env,db):

        # For each donor, the costs need to be input more specifically
        # before inputting into the different databases

        # CONSUMABLES

        cost_ET_isolation = sum(db.unit_planar_ets_cost*self.planar_ets_spent_isolation)

        cost_ET_expansion = sum(db.unit_planar_ets_cost*self.planar_ets_spent_expansion)

        cost_cryovials = self.cryovials_spent*db.cost_cryovial

        # REAGENTS

        cost_cm_isolation = self.volume_media_spent_isolation*db.cost_dmem_fbs

        cost_tryple_isolation = db.cost_tryple*self.volume_tryple_spent_isolation

        cost_fp_isolation = self.volume_FP_spent*db.cost_FP

        cost_pbs_isolation = self.volume_PBS_spent*db.cost_pbs

        cost_reagents_isolation = cost_cm_isolation + cost_tryple_isolation + cost_fp_isolation + cost_pbs_isolation

        cost_cm_expansion = self.volume_media_spent_expansion*db.cost_dmem_fbs

        cost_tryple_expansion = db.cost_tryple*self.volume_tryple_spent_expansion

        cost_reagents_expansion = cost_cm_expansion + cost_tryple_expansion

        cost_cm_dsp = self.volume_media_spent_DSP*db.cost_basal_medium

        cost_cryo_dsp = self.cryomedium_spent*db.cryomedium

        cost_reagents_DSP = self.volume_media_spent_DSP*db.cost_basal_medium + self.cryomedium_spent*db.cryomedium

        # DEPRECIATION AND FACILITY OPERATION COSTS

        cost_building_isolation = db.daily_facility_cost*math.ceil(self.time_isolation)

        cost_building_expansion = db.daily_facility_cost*math.ceil(self.time_expansion)

        cost_building_DSP = db.daily_facility_cost*math.ceil(self.time_DSP)

        cost_equipment_isolation = (db.daily_incubator_cost + db.daily_bsc_cost + db.daily_centrifuge_cost)*math.ceil(self.time_isolation)

        cost_equipment_expansion = (db.daily_incubator_cost + db.daily_bsc_cost + db.daily_centrifuge_cost)*math.ceil(self.time_expansion)

        cost_equipment_DSP = (db.daily_bsc_cost + db.daily_centrifuge_cost)*math.ceil(self.time_DSP)

        # LABOR

        cost_labor_isolation = (db.daily_worker*db.total_workers)*math.ceil(self.time_isolation)

        cost_labor_expansion = (db.daily_worker*db.total_workers)*math.ceil(self.time_expansion)

        cost_labor_DSP = (db.daily_worker*db.total_workers)*math.ceil(self.time_DSP)

        # RELEASE TESTING

        self.total_costs_release = db.cost_release_testing

        # TOTAL COSTS

        self.total_costs_isolation = cost_ET_isolation+cost_reagents_isolation

        self.total_costs_expansion = cost_ET_expansion+cost_reagents_expansion

        self.total_costs_DSP = cost_cryovials+cost_reagents_DSP

        self.total_costs_ets = cost_ET_expansion + cost_ET_isolation + cost_cryovials

        self.total_costs_reagents = cost_reagents_isolation + cost_reagents_expansion + cost_reagents_DSP

        self.total_hpl_costs = cost_cm_isolation + cost_cm_expansion

        self.total_costs_donor = self.total_costs_isolation + self.total_costs_expansion + self.total_costs_DSP + self.total_costs_release

        # print('Total costs of donor inside donor %.2f' % self.total_costs_donor)