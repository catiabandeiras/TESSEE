import numpy as np
import pandas as pd
from numpy.random import *

from individual import *

from databases_heor import *

number_runs = #INSERT NUMBER OF PATIENTS TO SIMULATE

#Initialize the structures to enable saving of results

years_simul = #INSERT YEARS OF SIMULATION

#INSERT THE NAMES OF HEALTH STATES AND CORRESPONDING INDEXES

dict_hs = {

}

#INSERT THE NAMES OF SPECIFIC HEALTH COMPLICATIONS AND CORRESPONDING INDEXES

dict_complic = {

}

#CREATE STANDARD OF CARE AND STEM CELL THERAPY RESULT DATABASES + INCREMENTAL COST EFFECTIVENESS RATIO (ICER)

soc_results = pd.DataFrame(data = np.zeros((number_runs,2)),
	index = range(number_runs),
	columns=['cost','qaly'])

stemcell_results = pd.DataFrame(data = np.zeros((number_runs,2)),
	index = range(number_runs),
	columns=['cost','qaly'])

icer = pd.DataFrame(data = np.zeros((number_runs,2)),
	index = range(number_runs),
	columns = ['icer','tag'])

icer_year = pd.DataFrame(data = np.zeros((number_runs,years_simul)),
	index = range(number_runs),
	columns = range(years_simul))


for run in range(number_runs):

	age = randint(18,35,size=1)

	indiv_soc = Individual('insulin',age[0],2,years_simul)

	##LEAVE IT HERE

	indiv_ins.run_indiv(years_simul)

	insulin_results.loc[run,'cost'] = round(sum(indiv_ins.costs_year),2)

	#print('qalys year insulin')
	#print(indiv_ins.qaly_year)

	insulin_results.loc[run,'qaly'] = round(indiv_ins.qalys,2)

	insulin_complic[run] = sum(indiv_ins.complications)

	insulin_complic_per_type.iloc[run,:] = indiv_ins.complications

	insulin_state_key.iloc[run,:] = indiv_ins.transplant_state_year

	#print('Beta cell analysis')

	indiv_sc = Individual('betacells',age[0],0,years_simul)
	
	indiv_sc.run_indiv(years_simul)

	transplant_results.loc[run,'cost'] = round(sum(indiv_sc.costs_year),2)



	
	#print('qalys year stem cells')
	#print(indiv_sc.qaly_year)

	transplant_results.loc[run,'qaly'] = round(indiv_sc.qalys,2)


	transplant_complic[run] = sum(indiv_sc.complications)

	transplant_complic_per_type.iloc[run,:] = indiv_sc.complications

	transplant_state.iloc[run,:] = indiv_sc.state_year

	transplant_state_key.iloc[run,:] = indiv_sc.transplant_state_year

	#print(len(np.arange(1,years_simul+1)))

	#if 

	#print((indiv_sc.utility_year - indiv_ins.utility_year)*np.arange(1,years_simul+1))
	#print(indiv_ins.costs_year)
	#print(indiv_ins.cum_costs_year)


	#print(indiv_sc.costs_year)
	#print(indiv_sc.cum_costs_year)

	#print(indiv_sc.costs_year - indiv_ins.costs_year)

	#print(indiv_sc.qaly_year - indiv_ins.qaly_year)

	if any(indiv_sc.qaly_year - indiv_ins.qaly_year == 0):

	#	print('Zero differences were found')

		icer_year.iloc[run,:] = (indiv_sc.cum_costs_year - indiv_ins.cum_costs_year)

	else:

		icer_year.iloc[run,:] = (indiv_sc.cum_costs_year - indiv_ins.cum_costs_year)/((indiv_sc.qaly_year - indiv_ins.qaly_year))

	#print(len(indiv_sc.costs_year))
	#print(len(indiv_sc.utility_year))

	#print(len(icer_year.iloc[run,:]))

#	icer_year.iloc[run,:] = (indiv_sc.costs_year - indiv_ins.costs_year)/((indiv_sc.utility_year - indiv_sc.utility_year))

#	if round(indiv_sc.qalys,2) <= round(indiv_ins.qalys,2):

#		break

	#print(indiv_sc.costs_year)

	#icer[run] = indiv.icer

# 	if round(indiv_sc.qalys,2) == round(indiv_ins.qalys,2):

# 		icer.iloc[run,0] = (round(sum(indiv_sc.costs_year),2) - round(sum(indiv_ins.costs_year),2))

# 		if icer.iloc[run,0] > 0:

# 			icer.iloc[run,1] = 'insulin'

# 		else:

# 			icer.iloc[run,1] = 'betacell'

# 	elif round(indiv_sc.qalys,2) > round(indiv_ins.qalys,2):

# 		icer.iloc[run,0] = (round(sum(indiv_sc.costs_year),2) - round(sum(indiv_ins.costs_year),2)) / (round(indiv_sc.qalys,2) - round(indiv_ins.qalys,2))

# 		if icer.iloc[run,0] > wtp_threshold:

# 			icer.iloc[run,1] = 'insulin'

# 		else:

# 			icer.iloc[run,1] = 'betacell'

# 	else:

# 		#Always assume more cost effectiveness from the insulin side

# 		icer.iloc[run,0] = (round(sum(indiv_sc.costs_year),2) - round(sum(indiv_ins.costs_year),2)) / (round(indiv_sc.qalys,2) - round(indiv_ins.qalys,2))

# 		icer.iloc[run,1] = 'insulin'

# 	#print('ICER of beta cell transplant is %d $/QALYs' % icer[run])

# print('ICER results')
# print(icer.describe())
# print(icer.quantile([0.025,0.975]))

qaly_difference = transplant_results.loc[:,'qaly'] - insulin_results.loc[:,'qaly']

#print(qaly_difference.describe())

#print(qaly_difference[qaly_difference < 0].describe())

insulin_results.to_csv('insulin_results.csv')

transplant_results.to_csv('transplant_results.csv')

insulin_complic.to_csv('insulin_complic.csv')

insulin_complic_per_type.to_csv('insulin_complic_per_type.csv')

insulin_state_key.to_csv('insulin_state_key.csv')

transplant_complic.to_csv('transplant_complic.csv')

transplant_state.to_csv('transplant_state.csv')

transplant_state_key.to_csv('transplant_state_key.csv')

transplant_complic_per_type.to_csv('transplant_complic_per_type.csv')

icer_year.to_csv('icer_year.csv')


print('Analysis 1 transplant, no immunosuppression, ends')

#icer.to_csv('icer.csv')

#print(insulin_results)
#print(transplant_results)