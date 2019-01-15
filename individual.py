import math
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from numpy.random import *
from database_diabetes import *

class Individual(object):

	def __init__(self,arm,init_age,init_state,years_simul):

		#Initialize the states to calculate the costs and QALYs

		#print(arm)

		self.age = init_age

		self.dead = 0

		self.initial_state = init_state

		self.current_state = init_state

		self.costs_year = np.zeros(years_simul)

		self.state_year = np.zeros(years_simul)

		self.transplant_state_year = np.zeros(years_simul)

		self.state_year[0] = init_state

		self.cum_costs_year = np.zeros(years_simul)

		self.qaly_year = np.zeros(years_simul)

		self.utility_year = np.zeros(years_simul)

		self.utility_year[0] = baseline_util[init_state]

		self.icer = 0

		self.life_years = 0

		self.qalys = 0

		self.complications = np.zeros(len(complic_probab_iit))

		self.complication_year = np.zeros(len(complic_probab_iit))

		self.transplant_complic = 1

		if arm == 'betacells':

			self.transplant_state = 0

		else:

			self.transplant_state = 2

		self.transplant_state_year[0] = self.transplant_state

	def run_indiv(self,years_simul):

		##Define the initial costs of the transplant

		##Assume that there is no initial graft failure


		if self.current_state < 2:

			#print('graft')

			self.costs_year[0] += cost_tr_proc + cost_device

			#Check the probability of initial complications 

			#is_init_compl = binomial(1,p = init_graft_complic_prob)

			is_init_compl = 0

			if is_init_compl == 1:

				self.current_state = 3

				self.transplant_complic = 1

				self.utility_year[0] += transplant_disut

		if self.current_state > 1:

			self.costs_year[0] += cost_iit_year

		#print(self.costs_year)

		self.cum_costs_year[0] = self.costs_year[0]

		self.qaly_year[0] = self.utility_year[0]

		years = 1

		while years < years_simul:

			while self.dead == 0:

				#Update the patient age

				self.update_state(years)

				if self.dead == 1:

					break

				else:

					years += 1

					self.age += 1

					if years == years_simul:

						break

			break

		# print('utilities tracking')
		# #print(years)
		# print(self.utility_year)
		# print(self.complications)
		# print(self.complication_year)
		# print(self.costs_year)

		# print('Final QALY sum')
		# print(sum(self.utility_year))

		self.life_years = years

		self.qalys = sum(self.utility_year)

		self.icer = sum(self.costs_year)/self.qalys

		# print('Additional life years %d ' % self.life_years)

		# print('Additional QALYs %.2f ' % self.qalys)

		# print('Final costs %d ' % sum(self.costs_year))

		# print('Final ICER is %d $/QALY' % self.icer)



	def update_state(self,years):

		#Update mortality state

		self.mortality()

		# print([years,self.dead])

		if self.dead == 0:

			self.graft(self.current_state,years)

			self.complication_prob_sampler(self.current_state,years)

			self.cost_per_state(self.current_state,years)

			self.utilities_per_state(self.current_state,years)

			self.state_year[years] = self.current_state

			self.transplant_state_year[years] = self.transplant_state

		else:

			self.transplant_state = 2

			self.costs_year[years] += 45000
			
			self.utility_year[years] = 0

			self.state_year[years:] = 4

			self.transplant_state_year[years:] = self.transplant_state

	def complication_prob_sampler(self,state,years):

		#Select the bimodality with the probabilities for each occurence

		#Get each state, irrespective of complications before or not

		if self.transplant_state == 0:

			vector_probs = complic_probab.iloc[0,:].values

		elif self.transplant_state == 1:

			vector_probs = complic_probab.iloc[1,:].values

		elif self.transplant_state == 2:

			vector_probs = complic_probab.iloc[2,:].values


		#vector_probs = [0.01,0.5]

		#print(vector_probs)
		ocurrences = binomial(1,vector_probs)

		is_complic = ocurrences.nonzero()[0]

		#Add to the number of complications not related to transplantation, but diabetes comorbidities

		#Limit the ocurrences to less than 2 in agreement with the Beckwith model. It was sampling too much

		if sum(self.complications) < 2 and sum(self.complications) + sum(ocurrences) <= 2:

			#Only add the complications once

			if is_complic.size > 0:

				list_of_complic = is_complic[0]

				if self.complications[list_of_complic] == 0:

					self.complications[list_of_complic] += ocurrences[list_of_complic]

					#Add the year it happened

					self.complication_year[ocurrences == 1] = years

		if sum(self.complications) > 0:

			self.current_state = 3 #State with complications


		# print('Complications sampling at year %d' % years)
		# print(self.complications)
		# print(self.complication_year)
		# print('The state after sampling is %d' % self.current_state)

	def cost_per_state(self, state, years):

		discount = discount_rate ** years

		if self.initial_state <2:

			self.costs_year[years] += cost_tr_year * discount

			#Inputting the costs of complications

			if self.transplant_complic == 1:

				self.costs_year[years] += cost_transplant_complic * discount

		if self.transplant_state == 2:

			#Inputting costs of IIT for people with dependence. In the end people with partial insulin independence, they have insulin independence but no levels of C-peptide

			self.costs_year[years] += cost_iit_year*discount

		#Then, look for the costs of complications (With or without transplant)

		if any(self.complications > 0):

			#Start with the costs for complications that happened the current year

			#print(baseline_costs[self.complication_year == years])

			self.costs_year[years] += sum(baseline_costs[self.complication_year == years])*discount

			#Add complications with a follow up time

			for comp in range(len(dict_diab_complic)):

				if self.complications[comp] > 0:

					difference = years - self.complication_year[comp]

					#print(difference)

					if difference < 5:

					#	print(follow_up_costs.iloc[int(difference)-1,comp])

						self.costs_year[years] += follow_up_costs.iloc[int(difference),comp]*discount

		#print(self.costs_year)

		self.cum_costs_year[years] = sum(self.costs_year[:years+1])

		#print(self.cum_costs_year)


	def utilities_per_state(self,state,years):

		#Make the decrement only associated with complications
		#We assume that the decrement applies only to the year where complications ocurred

		#Reduce the QALYs by discounting - doesn't need to multiply since each year deals with same discount rate!

		#discount = discount_rate ** years

		if any(self.complications > 0):

			self.utility_year[years] = (self.utility_year[years-1] + sum(complic_disutil[self.complication_year == years])) / discount_rate

		else:

			self.utility_year[years] = self.utility_year[years-1] / discount_rate



		self.qaly_year[years] = sum(self.utility_year[:years+1])


	def mortality(self):

		#Start with the baseline probability of death from the tables
		#Tables are now in the baseline, it's more correct

		prob_death = prob_death_db.iloc[self.age,1]

		#print('probability of death before complications %.5f' % prob_death)

		#Add the complications

		if any(self.complications[0:2] > 0):

			prob_death *= 2.4

		if any(self.complications[2:] > 0):

			prob_death *= 7.16

		#Update the state of life based on the probability

		#print('probability of death after complications assessment %.5f' % prob_death)

		self.dead = binomial(1,prob_death)

		#print('death %d' % self.dead)

	def graft(self,state,years):

		#Computes the expected behavior of graft failures

		#print(self.current_state)

		# if self.current_state == 0:

		# 	prob_fail = prob_graft_failure.iloc[years,:]

		# 	# print('Probability of remaining in state when full function')
		# 	# print(prob_fail)

		# 	state_graft = multinomial(1,pvals = prob_fail)

		# 	# print('State graft')
		# 	# print(state_graft)

		# 	if state_graft[1] == 1:

		# 		self.current_state = 1

		# 	elif state_graft[2] == 1:

		# 		self.current_state = 2		

		# if self.current_state == 1:

		# 	prob_fail = prob_graft_failure.iloc[years,1:]/sum(prob_graft_failure.iloc[years,1:])

		# 	# print('Probability of remaining in state when partial function')
		# 	# print(prob_fail)

		# 	state_graft = multinomial(1,pvals = prob_fail)

		# 	# print('State graft')
		# 	# print(state_graft)

		# 	if state_graft[-1] == 1:

		# 		self.current_state = 2

			# else:

			# 	#If the implant doesn't fail completely in the first year, it will achieve insulin independency
			# 	#Get back to becoming fully functional

			# 	if years == 1:

			# 		self.current_state = 2
			#print(self.current_state)

		#Start by the patients that do not have the failure

		if self.transplant_state != 2:

			#print(years)

			entry_state = self.current_state

			prob_fail = prob_full_graft_failure.iloc[years,1]

			is_full_fail = binomial(1,prob_fail)

			if is_full_fail == 1:

				self.transplant_state = 2

				if any(self.complications > 0):

					self.current_state = 3

				else:

					self.current_state = 2

			elif self.transplant_state == 0:

				prob_partial = prob_partial_graft_failure.iloc[years,1]

				is_partial = binomial(1,prob_partial)

				if is_partial == 1:

					self.transplant_state = 1

					if any(self.complications > 0):

						self.current_state = 3

					else:

						self.current_state = 1

				else:

					self.transplant_state = 0

					if any(self.complications > 0):

						self.current_state = 3

					else:

						self.current_state = 0

			elif self.transplant_state == 1:

				if any(self.complications > 0):

					self.current_state = 3

				else:

					self.current_state = 1


			# if entry_state == 0:

			# 	if state_graft[0] == 1:

			# 		self.current_state = 0

			# 	if state_graft[1] == 1:

			# 		self.current_state = 1

			# 	elif state_graft[2] == 1:

			# 		self.current_state = 2		

			# elif entry_state == 1:

			# 	if state_graft[2] == 0:

			# 		self.current_state = 1

			# 	else:

			# 		self.current_state = 2
