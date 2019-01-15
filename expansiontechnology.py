'''Expansiontechnology.py - File that includes the attributes associated to each individual expansion technology.
It can inherit attributes from the donor class.'''

import sys
import simpy
import math
#from tkinter import *
#from tkinter import ttk
import numpy as np
import pandas as pd
from numpy.random import *
from scipy.integrate import odeint

#from databases_esc import *

from databases_ipsc import *

class ExpansionTechnology(object):

	def __init__(self,name,sd,donor,env,db,number_wells):

		self.name = name

		if name == '6well':

			self.area = 19.6*6

		elif name == '500spinner':

			self.area = 500 #Just a dummy variable

		else:

			self.area = db.area_planar_ets[name]

		#print('Area of chosen expansion technology %d cm2' % self.area)

		self.seeding_density = sd

		#print('ET seeding density %d cells/cm2' % self.seeding_density)

		self.Xv0 = 0

		self.inc_time = donor.inc_time_passage_exp[donor.passage_no]

		self.time_to_confluence = donor.time_to_confluence[donor.passage_no]

		self.final_cells = 0

		self.expanded = 0

		#self.intra_donor_hd_factor = abs(normal(1,0.11)) #From Heathman for FBS
		#self.intra_donor_hd_factor = abs(normal(1,0.18)) #From Heathman for hPL

		#self.intra_donor_gr_factor = abs(normal(1,0.10)) #From Heathman for FBS
		#self.intra_donor_gr_factor = abs(normal(1,0.08)) #From Heathman for hPL

		self.intra_donor_hd_factor = 1
		self.intra_donor_gr_factor = 1

		self.growth_rate = donor.growth_rates[donor.passage_no]*self.intra_donor_gr_factor

		#print('Growth rate for expansion')
		#print(self.growth_rate)

		#self.harvest_dens = donor.base_hd[donor.passage_no]*self.intra_donor_hd_factor

		#Generate an harvest yield that can be modified with stochasticity. Use the triangular distribution such as Jenkins et al (2015)

		#self.harvest_yield = triangular(0.85,0.9,0.95)

		self.harvest_yield = 0.9

		#print('Growth rate of donor %d is %.2f at passage %d' % (donor.donor_index,self.growth_rate,donor.passage_no))
		#print('Harvest density of donor %d is %.2f at passage %d' % (donor.donor_index,self.harvest_dens,donor.passage_no))

		self.current_diff_stage = 0

		self.number_wells = number_wells

	def expansion(self,donor,env,db):

		#Expand the given ET type

		#Start by providing the initial cell number of this ET

		self.Xv0 = math.floor(self.seeding_density*db.area_planar_ets[self.name])

		#print('Seeding density of flask is %.2f' % self.seeding_density)
		#print('Seeded cells in a flask with %d initial cells at %.4f days' % (self.Xv0,env.now))

		#Seeding of the flask, get the volume media spent

		donor.volume_media_spent_expansion += db.media_volume_planar_ets[self.name]

		#print('Culture medium spent at time %d is %d' % (env.now,donor.volume_media_spent_expansion))

		#print('Media spent at expansion at time %d is %d' % (env.now,donor.volume_media_spent_expansion))

		yield env.timeout(db.et_operation_times_planar.loc[self.name,'seeding'])

#		while self.inc_time < self.time_to_confluence:

		while self.inc_time < donor.time_to_confluence[donor.passage_no]:
#		while self.final_cells < donor.base_hd[donor.passage_no]*db.area_planar_ets[self.name]/self.harvest_yield or self.final_cells < (donor.doses_per_donor*donor.cells_per_dose + db.cell_concentration*db.cryovial_vol)/len(donor.et_list):

			#Modify the harvesting density per flask with the intra donor variability

			#harvest_dens = donor.base_hd[donor.passage_no]*db.area_planar_ets[self.name]*self.intra_donor_hd_factor

			#print('Harvesting density %d' % harvest_dens)

			#if self.final_cells > harvest_dens:

			#	print('Something is wrong')

			#	break

			#Calculate the number of cells at this point

			yield env.timeout(db.time_bw_feedings)

			self.inc_time += db.time_bw_feedings

			#In

			if self.inc_time > donor.lag_phase_exp:

				self.final_cells = self.Xv0*math.exp(self.growth_rate*(self.inc_time-donor.lag_phase_exp))

			else:

				self.final_cells = self.Xv0

			#if self.inc_time > donor.time_to_confluence[donor.passage_no]: or self.final_cells > harvest_dens:

			# print('Incubation time %d days' % self.inc_time)
			# print('Cells at the time %d' % self.final_cells)
			# print('Current total elapsed time %d days' % env.now)

			if self.inc_time >= donor.time_to_confluence[donor.passage_no]:
				#self.final_cells = harvest_dens+1

			#print('Media spent at expansion at time %d is %d' % (env.now,donor.volume_media_spent_expansion))

			#print('Cells at time point %d are %d' % (env.now,self.final_cells))

			#yield env.timeout(db.et_operation_times_planar.loc[self.name,'feeding'])

			#if self.final_cells > harvest_dens:

				#print('Reached the number of cells required at %d' %env.now)

				break


			# if self.inc_time >= donor.time_to_confluence[donor.passage_no]:

			# 	print('Incubation time of the flask %d days' % self.inc_time)
			# 	print('Confluence time passed! Break culture for safety')
			# 	break

			else:

				yield env.timeout(db.et_operation_times_planar.loc[self.name,'feeding'])

				donor.volume_media_spent_expansion += db.media_volume_planar_ets[self.name]

				#print('Culture medium spent at time %d is %d' % (env.now,donor.volume_media_spent_expansion))

				self.inc_time += db.et_operation_times_planar.loc[self.name,'feeding']

			#print('Media spent at expansion at time %d after feeding is %d' % (env.now,donor.volume_media_spent_expansion))

			#print(self.cells_per_passage[self.passage_no])

		#Generates the harvesting part

		donor.volume_accutase_spent_expansion += db.harvesting_volume_planar_ets[self.name]

		#TrypLE requires a certain amount of media to inactivate it. Thus far, assume it is 2 times the volume of harvesting medium

		donor.volume_media_spent_expansion += 2*db.harvesting_volume_planar_ets[self.name]

		yield env.timeout(db.harvesting_time_incubator)

		#Add the cells to the total number of cells per passage

		#print(self.final_cells)
		#print(donor.base_hd[donor.passage_no]*db.area_planar_ets[self.name])

		#print('Final cells before flask harvesting are %d' % self.final_cells)

		cells_after_harvesting = self.final_cells*self.harvest_yield

		#print('Final cells after flask harvesting are %d' % cells_after_harvesting)

		donor.cells_per_passage_exp[donor.passage_no] += cells_after_harvesting

		self.expanded = 1

		donor.harvested_per_passage_exp[donor.passage_no] += 1

		#print('Culture media spent after harvesting at time %d is %d' % (env.now,donor.volume_media_spent_expansion))

	def growthcurve(self):

		#Growth curve model for the cell growth with the intra donor variability correction

		miu_max = self.growth_rate*self.intra_donor_gr_factor

		#use the time t as final for curve integration!

		t = np.linspace(0,self.inc_time,num=100)

		cellgrowth_int_handle = lambda Xv,t: self.growthcurve_ode(Xv,t,miu_max)

		cell_per_time = odeint(cellgrowth_int_handle,self.Xv0,t)

		#Returned cell numbers have to be integers

		return math.floor(cell_per_time[-1,-1])

	def growthcurve_ode(self,Xv,t,k):

		'''The most simple method for cell growth'''

		return k*Xv

	def differentiation(self,donor,env,db):

		#Determine what is the initial number of ESCs seeded in each 6 well plate

		if donor.diff_scheme == 'planar':

			volume_required = db.media_volume_well * self.number_wells

			volume_accutase = db.accutase_volume_well * self.number_wells

		elif donor.diff_scheme == 'spinner':

			volume_required = db.media_volume_spinner

			volume_accutase = db.accutase_volume_spinner

		#print('Start aggregation stage at %d days' % env.now)

		self.Xv0 = math.floor(self.seeding_density*volume_required)

		#Spend a certain amount of aggregation medium upon seeding

		#donor.volume_media_spent_aggregation += self.number_wells*db.media_volume_well

		#Start the feeding and then advance it all. Assume that the volume is always the same for every day that the medium is changed.

		donor.volume_media_spent_aggregation += volume_required*db.time_aggregation

		#print('Aggregation culture medium spent at time %d is %d' % (env.now,donor.volume_media_spent_aggregation))

		yield env.timeout(db.time_aggregation)

		# print('Finish aggregation stage at %d days' % env.now)

		# print('Start differentiation stage at %d days' % env.now)

		#Start the differentiation stage for good

		diff_time = 0

		self.current_diff_stage = 1

		final_diff_time = db.time_stage.values[-1]

		for stage in range(1,db.no_diff_stages + 1):

		#while diff_time < final_diff_time:

			#Select what is the number of the stage

			#self.diff_stage_selector(donor,env,db,diff_time)

			#According to the stage, select how many feeding schemes are conducted

			donor.volume_media_spent_differentiation[stage-1] += volume_required*db.number_feedings_stage[stage-1]

			#print('Differentiation culture medium of stage %d spent at time %d is %d' % (self.current_diff_stage,env.now,donor.volume_media_spent_differentiation[self.current_diff_stage-1]))

			yield env.timeout(db.time_stage[stage-1] - diff_time)

			#print('Current total simulation time %d days' % env.now)

			diff_time += (db.time_stage[stage-1]-diff_time)

			#print('Differentiation time %d days' % diff_time)

			#print('Stage %d ended' % stage)

			self.current_diff_stage += 1

			#if diff_time >= final_diff_time:

			#	break

		donor.volume_accutase_spent_differentiation += volume_accutase

		donor.harvested_diff += 1

		#print('Harvested differentiation flasks are %d' % donor.harvested_diff)

		donor.cells_per_passage_dif[donor.max_passages_dif-1] += self.Xv0 #Eliminate this and assume that the differentiation yield into cells of interest only comes at the DSP part
 
		self.differentiated = 1


	def diff_stage_selector(self,donor,env,db,diff_time):

		#Select what is the number of the current stage of differentiation

		number_stage = 0

		while number_stage < db.no_diff_stages:

			print('Current number of stage %d' % (number_stage+1))

			print('Current differentiation time %d' % diff_time)

			print('Time of differentiation stage %d' % db.time_stage.values[number_stage])

			if diff_time < db.time_stage.values[number_stage]:

				self.current_diff_stage = number_stage + 1

				print('Differentiation is in stage %d at %d days global, differentiation time %d days' % (self.current_diff_stage,env.now,diff_time))

				break

			else:

				number_stage += 1




