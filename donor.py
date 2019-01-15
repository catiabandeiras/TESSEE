'''donor.py - class responsible for the coordination of the donors'''

import sys
import simpy
import random
import math
import logging
import time
#from tkinter import *
#from tkinter import ttk
import numpy as np
import pandas as pd
from scipy.integrate import odeint
#Import what is going to solve the numpy validation
from numpy.random import *
from scipy import integrate, interpolate
from scipy import optimize

#Import the aditional modules

from databases import *

from expansiontechnology import *

class Donor(object):

	def __init__(self,env,donor_index,db):

		self.env = env

		self.donor_index = donor_index

		self.passage_no = 0

		self.isolated = 0

		self.expanded = 0

		self.finished = 0

		self.max_passages_exp = 5 #How to evaluate the maximum number of passages in pluripotent stem cells?

		self.max_passages_dif = 1 #Does it require passaging at any point?

		self.doses_per_donor = 1

		self.cells_per_dose = 75e6

		#Record the different cells of the passage

		self.cells_per_passage_exp = [0]*(self.max_passages_exp+1)
		self.cells_per_passage_dif = [0]*(self.max_passages_dif+1)

		#self.cells_per_passage[0] = self.initial_mnc

		#isolation_time = round(normal(9,1),0) #hPL
		#isolation_time = round(normal(12,1),0) #FBS

		#isolation_time = 12

		#For the Castiglia case, consider the P4 growth rate = 0 so that it is not considered
		#self.growth_rates = [0,0.073,0.259,0.3758,0]

		self.base_sd_exp = 50000
		self.base_sd_diff = 1e6

		#Initialize the simulations with the growth rate and base densities at the same as expansion

		self.base_sd = self.base_sd_exp


		self.growth_rate_exp = 0.3 #Start assuming a fixed growth rate for all the passages of expansion, add stochasticity later
		self.growth_rate_dif = 0.3 #Same here

		#self.growth_rates = [0,0.14,0.11,0.41,0]

		self.Xv0 = 0

		self.inc_time_passage_exp = [0]*(self.max_passages_exp+1)

		self.harvested_per_passage_exp = [0]*(self.max_passages_exp+1)

		self.et_list = []

		self.processed = 0

		self.cryovials_stored = 0

		self.tested = 0

		self.last_passage_name = ''

		self.last_passage_number = 0

		#db.occupied_inc += 1

		#Initialize media expenditures per donor

		self.volume_media_spent = 0

		self.volume_hr_spent = 0 # hr = harvesting reagent

		self.planar_ets_spent = pd.Series(data = [0]*len(db.names_of_planar_ets),
			index = db.names_of_planar_ets)


		self.volume_media_spent_expansion = 0

		self.volume_tryple_spent_expansion = 0

		self.planar_ets_spent_expansion = pd.Series(data = [0]*len(db.names_of_planar_ets),
			index = db.names_of_planar_ets)

		self.volume_media_spent_DSP = 0

		self.volume_PBS_spent = 0

		self.cryovials_spent = 0

		self.cryomedium_spent = 0

		self.initial_time_processing = 0

		#self.end_isolation_time = 0

		self.end_expansion_time = 0

		self.end_differentiation_time = 0

		self.end_DSP_time = 0

		self.time_expansion = 0

		self.time_differentiation = 0

		self.time_processing = 0

		self.time_DSP = 0

		self.total_costs_reprogramming = 0 #It's 0 for ESCs expansion, but some value when buying a vial of iPSCs

		self.total_costs_isolation = 0

		self.total_costs_expansion = 0

		self.total_costs_differentiation = 0

		self.total_costs_DSP = 0

		self.total_costs_release = 0

		self.total_costs_ets = 0

		self.total_costs_reagents = 0

		self.total_costs_donor = 0

		self.total_hpl_costs = 0

		#Input the costs of release testing as a single entity, since it is one dose per donor (autologous setting)

		# random_P0 = 0.50e6

		# while random_P0 <= 0.5e6:

		# 	random_P0 = normal(4.55e6,1.25e6)

		#random_P0 = 2e6

		#print(random_P0)

		initial_pscs = 1e6 #Make a dummy variable for either initial ESCs or the iPSCs coming from a vial

		self.cells_P0 = initial_pscs

		#self.cells_P0 = 982183

		#self.cpds = 0 #Are CPDs even relevant for iPSCs expansion?

		self.doses_no = 0

	def processing(self,env,db):

		'''Function that coordinates the complete processing of each donor'''

		self.initial_time_processing = env.now

		#Assuming that reprogramming is a service, in case it happens, launch the expansion process

		exp = self.expansion(env,db)

		env.process(exp)

		while self.expanded == 0:

			yield env.timeout(0.0001)

		# # # yield env.timeout(self.exp_time)

		print('Expansion of donor %d is finished at %d' % (self.donor_index,env.now))

		##Include the differentiation protocol

		dif = self.differentiation(env,db) #add the process

		env.process(dif)

		while self.differentiated == 0: #add the differentiation potential

			yield env.timeout(0.0001)

		# # # yield env.timeout(self.exp_time)

		print('Differentiation of donor %d is finished at %d' % (self.donor_index,env.now))

		# # #Launch the downstream process

		dsp = self.downstream(env,db)

		env.process(dsp)

		while self.finished == 0:

			yield env.timeout(0.0001)

		#yield env.timeout(self.dsp_time)

		print('Downstream processing of donor %d is finished after %d days' % (self.donor_index,env.now))

		qcr = self.release(env,db)

		env.process(qcr)

		while self.tested == 0:

			yield env.timeout(0.0001)

		print('Times spent in each unit operation:')
		print('Expansion: %d days' % self.time_expansion)
		print('Differentiation: %d days' % self.time_differentiation)
		print('DSP: %d days' % self.time_DSP)

		print('Processing of donor %d is finished after %d days' % (self.donor_index,env.now))

		self.time_processing = self.time_expansion + self.time_differentiation + self.time_DSP + db.qc_time

		self.processed = 1

		db.processed_donors += 1

		#db.occupied_inc -= 1

		db.total_doses += self.doses_per_donor

		#Sum all the expenses

		##About the expenses, they have to be checked beforehand, since the media for expansion and differentiation are not the same

		#db.volume_media_spent += self.volume_media_spent_expansion + self.volume_media_spent_DSP

		#db.volume_tryple_spent += self.volume_tryple_spent_isolation + self.volume_tryple_spent_expansion

		#self.planar_ets_spent = self.planar_ets_spent_isolation + self.planar_ets_spent_expansion

		#print(self.planar_ets_spent_isolation)

		#print(self.planar_ets_spent_expansion)

		#print(self.planar_ets_spent)

		#db.planar_ets_spent += self.planar_ets_spent_isolation + self.planar_ets_spent_expansion

		#db.volume_FP_spent += self.volume_FP_spent

		#db.cryovials_spent += self.cryovials_spent

		#db.cryomedium_spent += self.cryomedium_spent

		#self.cost_calculator_phase(env,db)

	def expansion(self,env,db):

		#Initialize the passage number and number of doses produced

		#passage_no = 1

		doses_no = 0

		# print('At P%d, the number of doses is %d' %(self.passage_no,doses_no))
		# print('At P%d, the number of cells is %d' %(self.passage_no,self.cells_per_passage[self.passage_no]))


		while self.passage_no <= self.max_passages and doses_no <= self.doses_per_donor:

			#Update the initial cells of the actual passage

			#self.cells_per_passage[self.passage_no] += self.cells_per_passage[self.passage_no-1]

			#Create the ETs

			self.create_ets(env,db,'expansion')

			for et in self.et_list:

				print('Launch')

				et_exp = et.expansion(self,env,db)

				env.process(et_exp)

				yield env.timeout(0)

			while self.harvested_per_passage[self.passage_no] < len(self.et_list):

				yield env.timeout(0.0001)

			#self.doses_no = math.floor(self.cells_per_passage[self.passage_no]/self.cells_per_dose)

			#print('Initial number of cells breakdown')
			#print(self.cells_P0)
			#print(self.initial_cells_for_expansion)

			self.cpds = round(math.log(self.cells_per_dose*self.doses_per_donor/self.initial_cells_for_expansion)/math.log(2),2)


			#self.cpds = round(math.log(self.cells_per_passage[self.passage_no]/self.initial_cells_for_expansion)/math.log(2),2)

			print('Final CPDs are %.2f' %self.cpds)

			# print('End of passage %d' %self.passage_no)
			# print('Number of cells %d' %self.cells_per_passage[self.passage_no])
			# print('Number of cells per dose %d' %self.cells_per_dose)
			# print('Number of doses %d' %doses_no)

			#print('Number of doses')
			#print(doses_no)

			#Account for loss of yield after downstream processing


			print('Comparisons before stopping condition')

			cells_for_yield = self.cells_per_passage[self.passage_no]
			#cells_for_stopping = self.doses_per_donor*self.cells_per_dose + db.cell_concentration*db.cryovial_vol
			cells_for_stopping = self.doses_per_donor*self.cells_per_dose

			print(cells_for_yield)
			print(cells_for_stopping)


			if cells_for_yield >= cells_for_stopping:

				#self.last_passage_name =

				print('Cells per passage reached at %d cells' % (self.cells_per_passage[self.passage_no]))
				print('Expected to be %d cells after DSP' % (self.cells_per_passage[self.passage_no]))

				break

			elif self.passage_no == self.max_passages:

				print('Could not reach the demand in time, but doses will be added with fewer cells')

				#self.doses_no = math.ceil(self.cells_per_passage[self.passage_no]/self.cells_per_dose)

				break

			else:

				self.passage_no += 1

		#Make the timeouts. So far, use only the fixed harvesting time in incubator.

		#print('Total number of cells of donor %d are %d' % (self.donor_index,self.cells_per_passage[self.passage_no]))

		#print('Expansion of donor %d finished at passage %d with %d doses produced' % (self.donor_index,self.passage_no,self.doses_no))

		self.expanded = 1

		self.end_expansion_time = env.now

		self.time_expansion = self.end_expansion_time - self.end_isolation_time

	def differentiation(self,env,db):

		#Method to go through the differentiation stages into pancreatic lineages.

		self.create_ets(env,db,'differentiation')

		initial_pscs = self.cells_per_passage[self.passage_no]

		differentiation_time = 0

		differentiation_time_end = db.time_stage[end]

		number_wells = initial_pscs/(db.media_volume_well * self.base_sd_diff)

		number_plates = number_wells / 6

		yield env.timeout(0.001)

	def downstream(self,env,db):

		#For simplicity, use the centrifuge as common just as a proof of concept
		#Use the fixed VR time

		print('Cells per passage before DSP are %d' % self.cells_per_passage[self.passage_no])

		washes_done = 0

		while washes_done < db.number_washes:

			#Needs to get the volume of the largest vessel used in the last passage
			#print('Parameters for DSP volume expenditures')
			#print(self.last_passage_name)
			#print(self.last_passage_number)

			self.volume_media_spent_DSP += self.last_passage_number*db.media_volume_planar_ets[self.last_passage_name] #Needs to be corrected!

			yield env.timeout(db.wash_vr_time)

			washes_done += 1

		#Input a final step yield for volume reduction and fill finish

		#self.cells_per_passage[self.passage_no] *= (db.vr_yield*db.ff_yield)

		#Put the cells in a culture medium to attain the right concentration

		print('Number of cells after DSP are %d' % self.cells_per_passage[self.passage_no])

		cryovials_no = math.ceil(self.cells_per_passage[self.passage_no]/(db.cell_concentration*db.cryovial_vol))

		print('Number of cryovials at DSP are %d' % cryovials_no)

		self.cryovials_spent += cryovials_no

		self.volume_media_spent_DSP += cryovials_no*db.cryovial_vol*(1-db.cryomedium_ratio)

		self.cryomedium_spent += cryovials_no*db.cryovial_vol*db.cryomedium_ratio

		self.cryovials_stored = cryovials_no

		yield env.timeout(db.finish_time)

		print('Cells processed from donor %d after DSP are %d' % (self.donor_index,self.cells_per_passage[self.passage_no]))

		self.end_DSP_time = env.now

		self.time_DSP = self.end_DSP_time - self.end_expansion_time

		self.finished = 1

	def release(self,env,db):

		'''Process related to passing the release or not'''

		#Sacrifice one of the vials for testing

		self.cryovials_stored -= 1

		yield env.timeout(db.qc_time)

		#Draw from a binomial distribution if it's a pass or a fail

		switch = np.random.binomial(1,db.pass_release_ratio)

		#print('qc test result')
		#print(switch)

		if switch == 1:

			print('Release testing of donor %d passed.' % self.donor_index)

			self.doses_no = math.floor(self.cryovials_stored * db.cell_concentration * db.cryovial_vol /self.cells_per_dose)

			#print(self.cryovials_stored)
			#print(db.cell_concentration)
			#print(db.cryovial_vol)

			print('Doses of donor %d are %d' % (self.donor_index,self.doses_no))

			self.doses_per_donor = self.doses_no

			print('Processing of donor %d finished at passage %d with %d doses produced' % (self.donor_index,self.passage_no,self.doses_no))

			print('Cells from donor %d are stored in %d cryovials with a concentration of %.2f cells/ml' % (self.donor_index,self.cryovials_stored,db.cell_concentration))

		else:

			print('Release testing of donor %d failed.' % self.donor_index)

			#Production of doses is cancelled. All goes to disposal!

			self.doses_per_donor = 0

		self.tested = 1

	def create_ets(self,env,db,stage):

		#First of all, choose what is the parametrization

		if stage == 'expansion':

			self.base_sd = self.base_sd_exp

			cells_to_seed = self.cells_per_passage[self.passage_no-1]

			print('Baseline seeding density %d' % self.base_sd)

			cells_per_et = self.base_sd*db.area_planar_ets

			no_ets = np.floor(cells_to_seed/cells_per_et)

			#sd_per_et = cells_to_seed/(db.area_planar_ets*no_ets)

			#Replace the infinite values by zero

			print('Number of expansion technologies generated')
			print(no_ets)

			#Choose the index

			et_name = ''
			sd_error = 1e16
			sd_new = 0

			#Choose, from options that yield the minimum expansion technologies, the ones that waste less cells

			#Subset first so that only values larger than zero appear. Sometimes the best option is not to have a zero

			no_ets_subset = no_ets[no_ets>0]
			#print(no_ets_subset)

			minimum_ets = min(no_ets_subset) #So that the minimum number of ETs is never zero!

			for name in db.names_of_planar_ets:

				#if no_ets[name] == min(no_ets):

				if no_ets[name] == minimum_ets:

					sd_per_et = cells_to_seed/(db.area_planar_ets[name]*minimum_ets)

					sd_error_aux = abs(sd_per_et-self.base_sd)

					if sd_error_aux < sd_error:

						sd_error = sd_error_aux

						et_name = name

						#sd_new = sd_per_et[name] #Old assumption

						sd_new = self.base_sd #Since no relationships between seeding density are incorporated, we will keep the seeding density

			#Create the list of expansion technologies

			print('name of the technology %s' % et_name)
			print(no_ets[et_name])

			#Make a cut on the number of cells

			#print('Base harvesting cells per flask')
			#print(self.base_hd[self.passage_no]*db.area_planar_ets[et_name])

			#Put at the minimum of harvesting yield! - 0.85! also account for the intra donor variability in density

			generated_cells = (self.base_hd[self.passage_no]*no_ets[et_name]*db.area_planar_ets[et_name])*0.85

			needed_cells = self.cells_per_dose*self.doses_per_donor + db.cell_concentration*db.cryovial_vol


			#if stage == 'expansion' and generated_cells > needed_cells:

				#generated_cells = self.base_hd[self.passage_no]*no_ets[et_name]*db.area_planar_ets[et_name]

				#needed_cells = self.cells_per_dose*self.doses_per_donor

				# print('The calculations are for %d cells and only %d are needed' % (generated_cells,needed_cells))

				# print('Number of %s before: %d' %(et_name,no_ets[et_name]))

				# no_ets[et_name] = math.ceil(needed_cells/(self.base_hd[self.passage_no]*db.area_planar_ets[et_name]/0.85))

				# print('Number of %s after: %d' %(et_name,no_ets[et_name]))

			if stage == 'expansion' and no_ets[et_name] > db.units_resource_planar.loc[et_name,'incubator']:

				no_ets[et_name] = db.units_resource_planar.loc[et_name,'incubator']

			elif generated_cells > needed_cells:

				no_ets[et_name] = math.ceil(needed_cells/(self.base_hd[self.passage_no]*db.area_planar_ets[et_name]*0.85)) #Contar com o downstream

			#sd_new = cells_to_seed/(no_ets[et_name]*db.area_planar_ets[et_name])

			self.et_list = [ExpansionTechnology(et_name,sd_new,self,env,db) for index in range(int(no_ets[et_name]))]

			#Add to the spent varieties of the donor

			self.planar_ets_spent_expansion[et_name] += int(no_ets[et_name])

			print('Passage %d of donor %d generates %d %s' % (self.passage_no,self.donor_index,no_ets[et_name],et_name))

			#Assume that this is the last passage until further information

			self.last_passage_name = et_name

			self.last_passage_number = no_ets[et_name]

			if self.passage_no == 1:

				self.initial_cells_for_expansion = round(no_ets[et_name]*db.area_planar_ets[et_name]*sd_new,0)

				print('Seeding density at P1 %d cells/cm^2' %sd_new)
				print('Initial cells for expansion seeded at P1 are %d' % self.initial_cells_for_expansion)

			# for et in self.et_list:

			# 	print('A %s was created with a seeding density of %d cells/cm^2 for passage %d' % (et.name,et.seeding_density,self.passage_no))

		elif stage == 'differentiation':

			self.base_sd = self.base_sd_diff

			cells_to_seed = self.cells_per_passage[self.passage_no]

			

	def cost_calculator_phase(self,env,db):

		cost_ET_isolation = sum(db.unit_planar_ets_cost*self.planar_ets_spent_isolation)

		#print('Costs of ETs in isolation: %.2f EUR' % cost_ET_isolation)

		#print('ETs spent in isolation')
		#print(self.planar_ets_spent_isolation)

		cost_ET_expansion = sum(db.unit_planar_ets_cost*self.planar_ets_spent_expansion)

		cost_cryovials = self.cryovials_spent*db.cost_cryovial

		cost_cm_isolation = self.volume_media_spent_isolation*db.cost_dmem_fbs

		cost_tryple_isolation = db.cost_tryple*self.volume_tryple_spent_isolation

		cost_fp_isolation = self.volume_FP_spent*db.cost_FP

		cost_pbs_isolation = self.volume_PBS_spent*db.cost_pbs

		cost_reagents_isolation = cost_cm_isolation + cost_tryple_isolation + cost_fp_isolation + cost_pbs_isolation

		#cost_cm_isol = self.volume_media_spent_isolation*db.cost_dmem_fbs

		print('Total culture media spent in isolation %d' % self.volume_media_spent_isolation)

		print('Total cost of culture media in isolation %.2f EUR' % cost_cm_isolation)

		print('Costs of reagents in isolation: %.2f EUR' % cost_reagents_isolation)

		cost_cm_expansion = self.volume_media_spent_expansion*db.cost_dmem_fbs

		cost_tryple_expansion = db.cost_tryple*self.volume_tryple_spent_expansion

		cost_reagents_expansion = cost_cm_expansion + cost_tryple_expansion

		print('Total culture media spent in expansion %d' % self.volume_media_spent_expansion)

		print('Total cost of culture media in expansion %.2f EUR' % cost_cm_expansion)

		print('Costs of reagents in expansion: %.2f EUR' % cost_reagents_expansion)

		print('ETs spent in expansion')
		print(self.planar_ets_spent_expansion)

		cost_cm_dsp = self.volume_media_spent_DSP*db.cost_basal_medium

		cost_cryo_dsp = self.cryomedium_spent*db.cryomedium

		cost_reagents_DSP = self.volume_media_spent_DSP*db.cost_basal_medium + self.cryomedium_spent*db.cryomedium

		print(cost_reagents_DSP)

		print('Total culture media spent in DSP %d' % self.volume_media_spent_DSP)
		print('Costs of reagents in downstream: %.2f EUR' % cost_reagents_DSP)

		cost_building_isolation = db.daily_facility_cost*math.ceil(self.time_isolation)

		cost_building_expansion = db.daily_facility_cost*math.ceil(self.time_expansion)

		cost_building_DSP = db.daily_facility_cost*math.ceil(self.time_DSP)

		#For the costs of equipment, we consider that only one piece of equipment is used in each donor expansion

		cost_equipment_isolation = (db.daily_incubator_cost + db.daily_bsc_cost + db.daily_centrifuge_cost)*math.ceil(self.time_isolation)

		cost_equipment_expansion = (db.daily_incubator_cost + db.daily_bsc_cost + db.daily_centrifuge_cost)*math.ceil(self.time_expansion)

		cost_equipment_DSP = (db.daily_bsc_cost + db.daily_centrifuge_cost)*math.ceil(self.time_DSP)

		#For the workers, do we consider that all the workers end up spending time in the tasks or not?

		cost_labor_isolation = (db.daily_worker*db.total_workers)*math.ceil(self.time_isolation)

		#print('Costs of labor in isolation: %.2f EUR' % cost_labor_isolation)

		#print('Costs of ETs in expansion: %.2f EUR' % cost_ET_expansion)

		#print('Costs of reagents in expansion: %.2f EUR' % cost_reagents_expansion)

		cost_labor_expansion = (db.daily_worker*db.total_workers)*math.ceil(self.time_expansion)

		#print('Costs of labor in expansion: %.2f EUR' % cost_labor_expansion)

		cost_labor_DSP = (db.daily_worker*db.total_workers)*math.ceil(self.time_DSP)

		#Input the costs of release testing as a single entity, since it is one dose per donor (autologous setting)

		self.total_costs_release = db.cost_release_testing

		#Takes off the total costs of building and equipment for now since they are variable. also take off the labor

		self.total_costs_isolation = cost_ET_isolation+cost_reagents_isolation

		self.total_costs_expansion = cost_ET_expansion+cost_reagents_expansion

		self.total_costs_DSP = cost_cryovials+cost_reagents_DSP

		self.total_costs_ets = cost_ET_expansion + cost_ET_isolation + cost_cryovials

		self.total_costs_reagents = cost_reagents_isolation + cost_reagents_expansion + cost_reagents_DSP

		self.total_hpl_costs = cost_cm_isolation + cost_cm_expansion

		self.total_costs_donor = self.total_costs_isolation + self.total_costs_expansion + self.total_costs_DSP + self.total_costs_release

		#self.total_costs_donor_dose = self.total_costs_donor/self.doses_per_donor

		#print('Total reagent costs %.2f' % self.total_costs_reagents)
		#print('Total hPL costs %.2f' % self.total_hpl_costs)


		#print('Total variable costs of donor %d are %d EUR' % (self.donor_index,self.total_costs_donor))
		#print('Total variable costs of donor %d per dose are %d EUR' % (self.donor_index,self.total_costs_donor_dose))
		#print('Costs of isolation: %.2f EUR'%self.total_costs_isolation)
		#print('Costs of expansion: %.2f EUR'%self.total_costs_expansion)
		#print('Costs of DSP: %.2f EUR'%self.total_costs_DSP)
		#print('Costs of release testing: %.2f EUR'%self.total_costs_release)


