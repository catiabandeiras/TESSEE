'''stemfactory.py - Class responsible for keeping the facility working'''

#Importing generic modules of Python
import sys
import simpy
import random
import math
import logging
import time
import pandas as pd

import matplotlib.pyplot as plt

#Call the other built-in modules - import the MACS module to start testing

from donor_esc_macs import *

#from databases_esc import *

from databases_ipsc import *

class StemFactory():

	def __init__(self,env,db):

		self.env = env

		self.cost_ET = 0

		self.cost_reagents = 0

		self.cost_building = 0

		self.cost_equipment = 0

		self.cost_labor = 0

		self.cost_qc = 0

		self.total_costs = 0

		self.costs_per_dose = 0

		self.total_doses = 0

		self.number_donors = 1 #Scale up process

		self.processed_donors = 0

		self.inc_list = simpy.Resource(env,db.total_incubators)

		self.bsc_list = simpy.Resource(env,db.total_bscs)

		self.cost_donor_stage = pd.DataFrame(
									data = np.zeros((self.number_donors,4)),
									index = range(self.number_donors),
									columns = ['expansion','differentiation','dsp','release']
									)

		self.times_process_stage = pd.DataFrame(
									data = np.zeros((self.number_donors,4)),
									index = range(self.number_donors),
									columns = ['expansion','differentiation','dsp','release']
									)

		self.cost_donor_variable = pd.DataFrame(
									data = np.zeros((self.number_donors,6)),
									index = range(self.number_donors),
									columns = ['ets','reagents','qc','building','equipment','labor']
									)

		self.time_processing_histogram = pd.Series(
			data = np.zeros((self.number_donors)),
			index = range(self.number_donors)
			)

		self.total_cells_dose = pd.Series(
			data = np.zeros((self.number_donors)),
			index = range(self.number_donors)
			)

		self.total_costs_dose = pd.Series(
			data = np.zeros((self.number_donors)),
			index = range(self.number_donors)
			)

		self.total_doses_per_donor = pd.Series(
			data = np.zeros((self.number_donors)),
			index = range(self.number_donors)
			)

		self.total_passages_donor = pd.Series(
			data = np.zeros((self.number_donors)),
			index = range(self.number_donors)
			)

		self.growth_rates_dataframe = pd.DataFrame(
									data = np.zeros((self.number_donors,6)),
									index = range(self.number_donors),
									columns = ['P0','P1','P2','P3','P4','P5']
									)

		# self.hd_dataframe = pd.DataFrame(
		# 							data = np.zeros((self.number_donors,6)),
		# 							index = range(self.number_donors),
		# 							columns = ['P0','P1','P2','P3','P4','P5']
		# 							)

		self.ets_per_donor = pd.DataFrame(
									data = np.zeros((self.number_donors,8)),
									index = range(self.number_donors),
									columns = db.names_of_planar_ets
									)

		self.P0_histogram = pd.Series(
			data = np.zeros((self.number_donors)),
			index = range(self.number_donors)
			)

		self.cpds_histogram = pd.Series(
			data = np.zeros((self.number_donors)),
			index = range(self.number_donors)
			)

		self.total_donor_costs_histogram = pd.Series(
			data = np.zeros((self.number_donors)),
			index = range(self.number_donors)
			)

		self.costs_failed_doses = 0


	def run(self,env,db):

		'''Define how the GMP facility operates'''

		#Initialize the donors

		#donor_index = 0

		#For autologous therapies, we just want to know how many donors we want to process
		#The annual demand compliance doesn't exist here

		#Create a list of donors

		#self.donors_list = [Donor(env,donor_index,db) for donor_index in range(self.number_donors)]

		for donor_index in range(self.number_donors):

			env.process(self.donor(env,donor_index,db))

		while db.processed_donors < self.number_donors:

			yield env.timeout(0.1)

		print('All donors processed after %.4f days!' % env.now)

		#Calculate the costs

		self.cost_calculator(env,db)

		#Generate the histogram

		#self.histogram_generator(env,db)
		# 		break



	def donor(self,env,donor_index,db):

		dn = Donor(env,donor_index,db)

		#print('Donor %d arrives to BSC at %.4f' % (donor_index,env.now))

		#with self.bsc_list.request() as request:

		#	yield request

		print('Donor %d arrives to incubator at %.4f' % (donor_index,env.now))

		with self.inc_list.request() as request:

			yield request

			print('Donor %d enters incubator at %.4f' % (donor_index,env.now))

			donor_proc = dn.processing(env,db)

			yield env.process(donor_proc)

			print('Donor %d leaves incubator at %.4f' % (donor_index,env.now))

			while dn.processed == 0:

				yield env.timeout(0.0001)

			#Sum the number of doses produced

			self.total_doses += dn.doses_no

			#Adds the variable costs of the donor to the dataframes

			self.cost_donor_stage.loc[dn.donor_index,'expansion'] = dn.total_costs_expansion
			self.cost_donor_stage.loc[dn.donor_index,'differentiation'] = dn.total_costs_differentiation
			self.cost_donor_stage.loc[dn.donor_index,'dsp'] = dn.total_costs_DSP
			self.cost_donor_stage.loc[dn.donor_index,'release'] = dn.total_costs_release

			self.times_process_stage.loc[dn.donor_index,'expansion'] = dn.time_expansion
			self.times_process_stage.loc[dn.donor_index,'differentiation'] = dn.time_differentiation
			self.times_process_stage.loc[dn.donor_index,'dsp'] = dn.time_DSP
			self.times_process_stage.loc[dn.donor_index,'release'] = db.qc_time

			self.cost_donor_stage.to_csv('cost_donor_stage.csv')

			# print('Check if things were well added')
			# print(self.cost_donor_stage)

			self.cost_donor_variable.loc[dn.donor_index,'ets'] = dn.total_costs_ets
			self.cost_donor_variable.loc[dn.donor_index,'reagents'] = dn.total_costs_reagents
			self.cost_donor_variable.loc[dn.donor_index,'qc'] = dn.total_costs_release

			self.growth_rates_dataframe.loc[dn.donor_index,:] = dn.growth_rates 

			#self.hd_dataframe.loc[dn.donor_index,:] = dn.base_hd

			self.cost_donor_variable.to_csv('cost_donor_variable.csv')

			self.P0_histogram[dn.donor_index] = dn.cells_per_passage_exp[0]/1e6

			self.cpds_histogram[dn.donor_index] = dn.cpds

			self.time_processing_histogram[dn.donor_index] = dn.time_processing

			self.total_donor_costs_histogram[dn.donor_index] = dn.total_costs_donor

			if dn.doses_no > 0:

				self.total_cells_dose[dn.donor_index] = dn.cells_after_dsp/(dn.doses_no*1e6)

				self.total_doses_per_donor[dn.donor_index] = dn.doses_no

			self.total_passages_donor[dn.donor_index] = dn.passage_no

			# print('On empty dataframe')
			# print(self.ets_per_donor)

			for et in db.names_of_planar_ets:
				self.ets_per_donor.loc[dn.donor_index,et] = dn.planar_ets_spent.loc[et]

			#self.ets_per_donor[dn.donor_index] = aux

			# print('On dataframe')
			# print(self.ets_per_donor)

			# self.ets_per_donor.to_csv('ets_per_donor.csv')

			# self.total_cells_dose.to_csv('total_cells_dose.csv')

			# self.total_donor_costs_histogram.to_csv('total_costs_donor.csv')

			# self.total_passages_donor.to_csv('total_passages_donor.csv')

			# P1_growthrates = self.growth_rates_dataframe.iloc[:,1]

			# print(P1_growthrates)

			# ax1 = P1_growthrates.plot.hist(bins=100)
			
			# plt.title('Distribution of P1 growth rates')
			# ax1.set_xlabel("Growth rates (day^-1)")
			# ax1.set_ylabel("Frequency")

			# fig1 = ax1.get_figure()
			# fig1.savefig('P1_growthrates.pdf')

			# self.growth_rates_dataframe.iloc[dn.donor_index,:] = dn.growth_rates

			# self.growth_rates_dataframe.to_csv('growth_rates.csv')

			# self.hd_dataframe.to_csv('harvesting_densities.csv')

			# self.cpds_histogram.to_csv('cpds.csv')

			# self.time_processing_histogram.to_csv('process_times.csv')

			# self.hPL_histogram.to_csv('hpl_costs.csv')


			# print('Check if things were well added')
			# print(self.cost_donor_variable)



				#while donor.processed == 0:

					#yield env.timeout(0.0001)

			#donor_index += 1

			#continue

			#while True:

		#If the cycle is finished and there are no more patients to treat, done

		# while True:

		# 	if db.processed_donors == self.number_donors:

		# 		print('All donors processed after %.4f days!' % env.now)

		# 		self.cost_calculator(env,db)

		# 		break

		# 	else:

		# 		yield env.timeout(0.1)

		# 		continue

	def cost_calculator(self,env,db):

		#print(db.unit_planar_ets_cost)

		#print(db.planar_ets_spent)

		#First, add all the fixed costs that vary with the total number of days in culture

		total_process_time = sum(self.time_processing_histogram)

		#print(total_process_time)

		self.cost_building = db.daily_facility_cost*math.ceil(env.now)

		self.cost_equipment = (db.daily_incubator_cost*db.total_incubators + db.daily_bsc_cost*db.total_bscs + db.daily_centrifuge_cost*db.total_centrifuges)*math.ceil(env.now)

		self.cost_labor = (db.daily_worker*db.total_workers)*math.ceil(env.now)

		print('Cost of building %d EUR' % self.cost_building)
		print('Cost of equipment %d EUR' % self.cost_equipment)
		print('Cost of labor %d EUR' % self.cost_labor)

		#Divide each of the costs by the total number of doses produced

		#Sum these costs to the total costs of each donor

		for i in range(self.number_donors):

			print('Process time of current donor')
			print(self.time_processing_histogram[i])
			#print(self.time_processing_histogram.iloc[i,1])

			self.cost_building_per_donor = self.cost_building*(self.time_processing_histogram[i]/total_process_time)

			#print('Cost of building per dose %d EUR' % self.cost_building_per_dose)

			self.cost_equipment_per_donor = self.cost_equipment*(self.time_processing_histogram[i]/total_process_time)

			self.cost_labor_per_donor = self.cost_labor*(self.time_processing_histogram[i]/total_process_time)

			self.cost_donor_variable.loc[i,'building'] = self.cost_building_per_donor
			self.cost_donor_variable.loc[i,'equipment'] = self.cost_equipment_per_donor
			self.cost_donor_variable.loc[i,'labor'] = self.cost_labor_per_donor

			self.total_donor_costs_histogram[i] += (self.cost_building_per_donor + self.cost_equipment_per_donor + self.cost_labor_per_donor)

			if self.total_doses_per_donor[i] == 0:

				self.costs_failed_doses += self.total_donor_costs_histogram[i]

			else:

				#If allogeneic!

				self.total_costs_dose[i] = self.total_donor_costs_histogram[i]/self.total_doses_per_donor[i]

				#Autologous

			
				self.total_costs_dose[i] = self.total_donor_costs_histogram[i]

			self.cost_donor_stage.loc[i,'expansion'] += (self.cost_building_per_donor+self.cost_equipment_per_donor+self.cost_labor_per_donor)*round(self.times_process_stage.loc[i,'expansion']/self.time_processing_histogram[i],2)
			self.cost_donor_stage.loc[i,'differentiation'] += (self.cost_building_per_donor+self.cost_equipment_per_donor+self.cost_labor_per_donor)*round(self.times_process_stage.loc[i,'differentiation']/self.time_processing_histogram[i],2)
			self.cost_donor_stage.loc[i,'dsp'] += (self.cost_building_per_donor+self.cost_equipment_per_donor+self.cost_labor_per_donor)*round(self.times_process_stage.loc[i,'dsp']/self.time_processing_histogram[i],2)
			self.cost_donor_stage.loc[i,'release'] += (self.cost_building_per_donor+self.cost_equipment_per_donor+self.cost_labor_per_donor)*round(self.times_process_stage.loc[i,'release']/self.time_processing_histogram[i],2)

			#print()
			
			#print(self.cost_donor_variable)
			#print(self.cost_building_per_dose)

			#print('Cost of building per dose %d EUR' % self.cost_donor_variable.loc[i,'building'])

		# print(self.times_process_stage)
		# print(self.times_process_stage.sum(1))
		# print(self.time_processing_histogram)

		#Add the costs of the failed doses equally distributed for all

		for i in range(self.number_donors):

			if self.total_doses_per_donor[i] > 0:

				#Allogeneic

				#self.total_costs_dose[i] += self.costs_failed_doses/sum(self.total_doses_per_donor)

				#Autologous

				self.total_costs_dose[i] += self.costs_failed_doses/self.number_donors

		self.time_processing_histogram.to_csv('process_times.csv')

		self.cost_donor_variable.to_csv('cost_donor_variable.csv')

		self.total_donor_costs_histogram.to_csv('total_costs_donor.csv')

		self.total_costs_dose.to_csv('total_costs_dose.csv')

		self.cost_donor_stage.to_csv('cost_donor_stage.csv')

		# print('Costs summed over variables')
		# print(self.cost_donor_variable.sum(1))

		# print('Costs summed over stages')
		# print(self.cost_donor_stage.sum(1))

		# print('Real total costs')
		# print(self.total_donor_costs_histogram)




		self.total_costs = sum(self.total_donor_costs_histogram)

		#self.costs_per_dose = sum(self.total_donor_costs_histogram)/sum(self.total_doses_per_donor)

		self.costs_per_dose = sum(self.total_donor_costs_histogram)

		print('Total costs of ETs %.2f EUR' %self.cost_ET)

		print('Total costs of reagents %.2f EUR' %self.cost_reagents)

		print('Total costs of building %.2f EUR' %self.cost_building)

		print('Total costs of equipment %.2f EUR' %self.cost_equipment)

		print('Total costs of labor %.2f EUR' %self.cost_labor)

		print('Total costs of quality control %.2f EUR' %self.cost_qc)

		print('Total costs of process %.2f EUR' % self.total_costs)

		print('Total costs per dose: %.2f EUR' % self.costs_per_dose)

	# def histogram_generator(self,env,db):

	# 	#Method used to generate the histograms for the distributions
	# 	#title('Number of P0 cells per donor')
	# 	# ax = self.P0_histogram.hist()
	# 	# fig = ax.get_figure()
	# 	# fig.savefig('hist_P0.pdf')
	# 	# plt.figure()
	# 	# self.P0_histogram.plot.hist()

	# 	plt.figure()

	# 	ax = self.P0_histogram.hist(bins=100)

	# 	#plt.xlim(0.5,10)
	# 	plt.title('Distribution of initial # of P0 BM-MSCs')
	# 	ax.set_xlabel("Million P0 BM-MSCs")
	# 	ax.set_ylabel("Frequency")

	# 	fig = ax.get_figure()
	# 	fig.savefig('hist_P0.pdf')

	# 	plt.figure()

	# 	self.P0_histogram.to_csv('P0_cells.csv')

	# 	ax = self.cpds_histogram.hist(bins=100)

	# 	plt.xlim(0,20)
	# 	plt.title('Distribution of CPDs')
	# 	ax.set_xlabel("CPDs")
	# 	ax.set_ylabel("Frequency")

	# 	fig = ax.get_figure()
	# 	fig.savefig('CPDs.pdf')

	# 	plt.figure()

	# 	self.cpds_histogram.to_csv('cpds.csv')

	# 	ax = self.time_processing_histogram.hist(bins=100)

	# 	#plt.xlim(0,100)
	# 	plt.title('Distribution of processing times')
	# 	ax.set_xlabel("Process times (days)")
	# 	ax.set_ylabel("Frequency")

	# 	fig = ax.get_figure()
	# 	fig.savefig('process_times.pdf')

	# 	plt.figure()

	# 	ax = self.total_donor_costs_histogram.hist(bins=100)

	# 	#plt.xlim(0,100)
	# 	plt.title('Distribution of total costs manufacturing/donor')
	# 	ax.set_xlabel("Processing costs (€)")
	# 	ax.set_ylabel("Frequency")

	# 	fig = ax.get_figure()
	# 	fig.savefig('total_donor_process_costs.pdf')

	# #def donor_launch(self,env,donor_index,db):

	# 	plt.figure()

	# 	ax = self.total_cells_dose.hist(bins=100)

	# 	#plt.xlim(0.5,10)
	# 	plt.title('Distribution of total cells/dose')
	# 	ax.set_xlabel("Million BM-MSCs")
	# 	ax.set_ylabel("Frequency")

	# 	fig = ax.get_figure()
	# 	fig.savefig('hist_final_cells.pdf')

	# 	plt.figure()

	# 	ax = self.cost_donor_stage.plot.hist(stacked=True, bins=100)
		
	# 	plt.title('Distribution of total costs manufacturing per stage')
	# 	ax.set_xlabel("Processing costs (€)")
	# 	ax.set_ylabel("Frequency")

	# 	fig = ax.get_figure()
	# 	fig.savefig('donor_costs_per_stage.pdf')

	# 	#def donor_launch(self,env,donor_index,db):

	# 	plt.figure()

	# 	ax = self.cost_donor_variable.plot.hist(stacked=True, bins=100)
		
	# 	plt.title('Distribution of total costs manufacturing per resource')
	# 	ax.set_xlabel("Processing costs (€)")
	# 	ax.set_ylabel("Frequency")

	# 	fig = ax.get_figure()
	# 	fig.savefig('donor_costs_per_resource.pdf')

		#Attempt to check what are the correspondences for the wrong donor costs
		#self.total_donor_costs_histogram.to_csv('total_costs_donor.csv')

		# P1_growthrates = self.growth_rates_dataframe.iloc[:,1]

		# print(P1_growthrates)

		# ax1 = P1_growthrates.plot.hist(bins=100)
		
		# plt.title('Distribution of P1 growth rates')
		# ax1.set_xlabel("Growth rates (day^-1)")
		# ax1.set_ylabel("Frequency")

		# fig1 = ax1.get_figure()
		# fig1.savefig('P1_growthrates.pdf')

		#self.growth_rates_dataframe.to_csv('growth_rates.csv')

