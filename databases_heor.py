import math
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from numpy.random import *

##MORTALITY RATES

prob_death_db = pd.read_csv('mortality_tables.csv',header=None)

prob_full_graft_failure = pd.read_csv('prob_full_graft_failure.csv',header=None)

prob_partial_graft_failure = pd.read_csv('prob_partial_graft_failure.csv',header=None)

#prob_full_graft_failure.iloc[1:,1] *= 0.25

#prob_partial_graft_failure.iloc[1:,1] *= 0.25

#print(prob_partial_graft_failure.head())

#Save all the costs, utilities in pandas databases

##INDEX AND COLUMNS DEFINITION

#Create a dictionary for the indexes

dict_diab_hs = {
	"full graft": 0,
	"partial graft": 1,
	"iit": 2,
	"complications": 3,
	"death": 4
}

#Create a dictionary for the complications

dict_diab_complic = {
	"hypoglicemia, in": 0,
	"hypoglicemia, er": 1,
	"stroke": 2,
	"pvd": 3,
	"ami": 4,
	"chf": 5,
	"amputation": 6,
	"ulcer": 7,
	"gangrene": 8,
	"esrd": 9,
	"blindness": 10,
	"edema": 11,
	"proliferative retinopathy": 12,
	"nonproliferative retinopathy": 13,
	"cataract": 14
}

##COMPLICATION DEFINITION

#Create a Pandas database for the probabilities of complications per disease

complic_probab = pd.DataFrame(data = np.zeros((len(dict_diab_hs),len(dict_diab_complic))),
	index = dict_diab_hs.keys(),
	columns = dict_diab_complic.keys()
	)

prob_hypo_in = 13.5
prob_hypo_er = 32.8

prob_stroke = 10.6
prob_pvd = 16.8
prob_ami = 9.1
prob_chf = 308

prob_amp = 5.7
prob_ulcer = 193.3
prob_gang = 250

prob_esrd = 19.2
prob_blind = 110
prob_edema = 110
prob_prol_ret = 50
prob_np_ret = 250
prob_cataract = 124

complic_probab_iit = np.zeros(len(dict_diab_complic))

complic_probab_iit[0] = prob_hypo_in
complic_probab_iit[1] = prob_hypo_er

complic_probab_iit[2] =prob_stroke
complic_probab_iit[3] =prob_pvd
complic_probab_iit[4] =prob_ami
complic_probab_iit[5] =prob_chf

complic_probab_iit[6] =prob_amp
complic_probab_iit[7] =prob_ulcer
complic_probab_iit[8] =prob_gang

complic_probab_iit[9] =prob_esrd
complic_probab_iit[10] =prob_blind
complic_probab_iit[11] =prob_edema
complic_probab_iit[12] =prob_prol_ret
complic_probab_iit[13] =prob_np_ret
complic_probab_iit[14] =prob_cataract

#Comes from incidences

complic_probab_iit = complic_probab_iit/10000

complic_probab_fg = complic_probab_iit*0.25
complic_probab_pg = complic_probab_iit*0.455

complic_probab.loc["full graft",:] = complic_probab_fg
complic_probab.loc["partial graft",:] = complic_probab_pg
complic_probab.loc["iit",:] = complic_probab_iit
complic_probab.loc["complications",:] = complic_probab_iit

init_graft_complic_prob = 0.65*0.6 #Based on Wallner et al, 2018

#print(complic_probab)

##UTILITIES DEFINITION



#Disutilities per complication
complic_disutil = pd.Series(data = np.zeros(len(dict_diab_complic)),
	index = dict_diab_complic.keys()
	)

disut_hypo_in = -0.047
disut_hypo_er = -0.004

disut_stroke = -0.291
disut_pvd = -0.11
disut_ami = -0.06
disut_chf = -0.058

disut_amp = -0.116
disut_ulcer = -0.083
disut_gang = 0

disut_esrd = -0.082

disut_blind = -0.024
disut_edema = -0.063
disut_prol_ret = -0.048
disut_np_ret = -0.048
disut_cataract = -0.063

complic_disutil[0] = disut_hypo_in
complic_disutil[1] = disut_hypo_er

complic_disutil[2] =disut_stroke
complic_disutil[3] =disut_pvd
complic_disutil[4] =disut_ami
complic_disutil[5] =disut_chf

complic_disutil[6] =disut_amp
complic_disutil[7] =disut_ulcer
complic_disutil[8] =disut_gang

complic_disutil[9] =disut_esrd
complic_disutil[10] =disut_blind
complic_disutil[11] =disut_edema
complic_disutil[12] =disut_prol_ret
complic_disutil[13] =disut_np_ret
complic_disutil[14] =disut_cataract

transplant_disut = -0.05*0.6

#Baseline utilities of each health state

baseline_util = pd.Series(data = np.zeros(len(dict_diab_hs)),
	index = dict_diab_hs.keys())

util_fg = 0.91
util_pg = 0.81
util_iit = 0.71

baseline_util[0] = util_fg
baseline_util[1] = util_pg
baseline_util[2] = util_iit

##COST DEFINITION

#Start with baseline new event costs

baseline_costs = pd.Series(data = np.zeros(len(dict_diab_complic)),
	index = dict_diab_complic.keys()
	)

cost_hypo_in = 19302
cost_hypo_er = 1524

cost_stroke = 13682
cost_pvd = 22527
cost_ami = 30181
cost_chf = 12952

cost_amp = 23825
cost_ulcer = 0
cost_gang = 0

cost_esrd = 0

cost_blind = 0
cost_edema = 0
cost_prol_ret = 0
cost_np_ret = 0
cost_cataract = 0

baseline_costs[0] = cost_hypo_in
baseline_costs[1] = cost_hypo_er

baseline_costs[2] =cost_stroke
baseline_costs[3] =cost_pvd
baseline_costs[4] =cost_ami
baseline_costs[5] =cost_chf

baseline_costs[6] =cost_amp
baseline_costs[7] =cost_ulcer
baseline_costs[8] =cost_gang

baseline_costs[9] =cost_esrd
baseline_costs[10] =cost_blind
baseline_costs[11] =cost_edema
baseline_costs[12] =cost_prol_ret
baseline_costs[13] =cost_np_ret
baseline_costs[14] =cost_cataract

baseline_costs *= 1.03 #2015 to 2017 USD

#Get the follow up costs per adverse event

follow_up_costs = pd.DataFrame(data = np.zeros((5,len(dict_diab_complic))),
	index = range(1,6),
	columns = dict_diab_complic.keys()
	)

#Cumulative follow up costs

follow_up_hypo_in = [184]*5
follow_up_hypo_ed = [185]*5

follow_up_stroke = [8719,6227,883,883,883]
follow_up_pvd = [9897,1466,368,368,368]
follow_up_ami = [14208,8633,81,81,81]
follow_up_chf = [10809,6041,3601,3601,3601]

follow_up_amp = [17371,1480,1480,1480,1480]
follow_up_ulcer = [2185,0,0,0,0]
follow_up_gang = [11378,1191,1191,1191,1191]

follow_up_esrd = [220187,220187,220187,220187,220187]

follow_up_blind = [2913,2913,2913,2913,2913]
follow_up_edema = [812,71,71,71,71]
follow_up_prol_ret = [625,71,71,71,71]
follow_up_np_ret = [630,256,256,256,256]
follow_up_cataract = [629,68,68,68,68]

follow_up_costs.iloc[:,0] = follow_up_hypo_in
follow_up_costs.iloc[:,1] = follow_up_hypo_ed

follow_up_costs.iloc[:,2] = follow_up_stroke
follow_up_costs.iloc[:,3] = follow_up_pvd
follow_up_costs.iloc[:,4] = follow_up_ami
follow_up_costs.iloc[:,5] = follow_up_chf

follow_up_costs.iloc[:,6] = follow_up_amp
follow_up_costs.iloc[:,7] = follow_up_ulcer
follow_up_costs.iloc[:,8] = follow_up_gang 

follow_up_costs.iloc[:,9] = follow_up_esrd

follow_up_costs.iloc[:,10] = follow_up_blind
follow_up_costs.iloc[:,11] = follow_up_edema
follow_up_costs.iloc[:,12] = follow_up_prol_ret
follow_up_costs.iloc[:,13] = follow_up_np_ret
follow_up_costs.iloc[:,14] = follow_up_cataract

follow_up_costs *= 1.03 #2015 to 2017 USD

##TRANSPLANTATION COSTS - Moassesfar 2016

cost_tr_proc = 23726.47*1.0676 #Inflon rate from 2012 to 2017

#cost_device = 100000/0.3

#cost_device = 100000

#Start running the analysis with the same costs of device as a transplant

#cost_device = 75468*1.0676

cost_device = 62470



cost_tr_year = 19000*1.181*0.6

cost_transplant_complic = 4900*1.181*0.6

##INSULIN ADMINISTRATION COSTS

cost_iit_year = 9601

##BASELINE DISCOUNT RATES

discount_rate = 1.03








