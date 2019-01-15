'''databases.py - contain all the classes associated with fixed values'''

import pandas as pd
import numpy as np
from numpy.random import *
#import simpy

class FacilityDB(object):

	'''Class containing all the associated variables to the construction of the stem factory'''

	def __init__(self):

		#Give the area in sq mt (can add later option for sq ft)

		self.facility_area = 400

		self.clean_room_ratio = 0.2

		#Costs in EUR for building the clean room and others

		self.cost_sqmt_clean_room = 5815

		self.cost_sqmt_other_room = 3392

		self.total_facility_cost = (self.cost_sqmt_clean_room*self.clean_room_ratio + self.cost_sqmt_other_room*(1-self.clean_room_ratio))*self.facility_area

		print('Total cost of facility is %.2f EUR' % self.total_facility_cost)

		#Facility depreciation period in days

		self.facility_depreciation_period = 15*365.25

		#Input the other costs related with gas supplies and other facilities

		co2_daily_cost = 6000/365.25

		gases_daily_cost = 15600/365.25

		add_supplies_daily_cost = 7900/365.25

		requal_daily_cost = 65400/365.25

		maint_daily_cost = 52800/365.25

		cleaning_daily_cost = 28000/365.25

		garment_daily_cost = 2000/365.25

		compilation_daily_cost = co2_daily_cost + gases_daily_cost + add_supplies_daily_cost + requal_daily_cost + maint_daily_cost + cleaning_daily_cost + garment_daily_cost 

		self.daily_facility_cost = self.total_facility_cost/self.facility_depreciation_period + compilation_daily_cost

		print('Daily cost of facility %.2f EUR' % self.daily_facility_cost)

		#Provide the resources
		#For autologous therapies, no 3D cultures - the goal is scale out

		self.total_workers = 3

		self.total_bscs = 4

		self.total_incubators = 8

		self.total_centrifuges = 4

		# self.total_workers = 9

		# self.total_bscs = 1

		# self.total_incubators = 1

		# self.total_centrifuges = 1

		self.occupied_inc = 0

		#Provide the costs of these resources

		self.daily_worker = 100

		self.unit_incubator_cost = 16000

		self.unit_bsc_cost = 15300

		self.unit_centrifuge_cost = 12000 #assume this from ThermoFisher

		#Equipment depreciation in days

		self.equipment_depreciation_period = 5*365.25

		self.daily_incubator_cost = self.unit_incubator_cost*self.total_incubators/self.equipment_depreciation_period

		self.daily_bsc_cost = self.unit_bsc_cost*self.total_bscs/self.equipment_depreciation_period

		self.daily_centrifuge_cost = self.unit_centrifuge_cost*self.total_centrifuges/self.equipment_depreciation_period

		#Input additional costs later, these suffice so far

		#Input the drivers associated with planar technologies, operating times, costs, etc

		self.names_of_planar_ets = ['t-flask25',
					   't-flask75',
					   't-flask175',
					   't-flask225',
					   'cellstack1',
					   'cellstack2',
					   'cellstack5',
					   'cellstack10']

		self.area_planar_ets = pd.Series(data = [25,75,175,225,636,636*2,636*5,636*10],
								index = self.names_of_planar_ets)

		self.media_volume_planar_ets = pd.Series(data = [5,15,35,45,130,260,650,1300],
									index = self.names_of_planar_ets)

		self.harvesting_volume_planar_ets = pd.Series(data = [1.75,3.5,7,9,25,50,125,250],
											index = self.names_of_planar_ets)

		self.units_resource_planar = pd.DataFrame(
									data = [[100,1],
											[100,1],
											[100,1],
											[100,1],
											[60,1],
											[60,1],
											[24,1],
											[12,1]],
									index = self.names_of_planar_ets,
									columns = ['incubator','worker']
									)

		self.unit_planar_ets_cost = pd.Series(data = [1.65,4.85,7.38,8.55,33.89,60.12,131.96,142.72],
												index = self.names_of_planar_ets)

		self.et_operation_times_planar = pd.DataFrame(
		    data = [[0.38/(24*10),0.38/(24*10),0.5/(24*10)],
		           [0.38/(24*10),0.38/(24*10),0.5/(24*10)],
		           [0.38/(24*10),0.38/(24*10),0.5/(24*10)],
		           [0.38/(24*10),0.38/(24*10),0.5/(24*10)],
		           [0.15/24,0.15/24,0.06/24],
		           [0.15/24,0.15/24,0.06/24],
		           [0.20/24,0.20/24,0.16/24],
		           [0.25/24,0.25/24,0.25/24]],
		    index = self.names_of_planar_ets,
		    columns = ['seeding','feeding','harvesting']
		    )

		self.harvesting_time_incubator = 14/(60*24)

		self.time_bw_feedings = 3

		#Get the costs of the reagents that are used in every step.
		#Start with DMEM + 10% FBS, then can work out for other materials

		#self.cost_dmem_fbs = 0.04*0.90 + 3.87*0.10 #hPL - from Heathman et al (2016) study, it's 10% hPLs

		self.cost_dmem_fbs = 0.04*0.90 + 2.12*0.10 #FBS

		self.cost_basal_medium = 0.04

		self.cost_tryple = 0.21

		#Initialize the volumes spent of each reagent

		self.volume_media_spent = 0

		self.volume_basal_medium_spent = 0

		self.volume_tryple_spent = 0

		self.planar_ets_spent = pd.Series(data = [0]*len(self.names_of_planar_ets),
											index = self.names_of_planar_ets)

		self.total_doses = 0

		self.processed_donors = 0

class IsolationDB(object):

	'''Class containing all the donor MSC isolation parameters'''

	def __init__(self):

		#self.isol_time = 4

		#Time to confluence to obtain the P0 cells

		#self.isol_conf_time = 14

		#The volume of BM received

		self.volume_BM = 10

		#Volume of Ficoll-Paque per tube (let's assume 50 ml)

		self.volume_FP = self.volume_BM*2 #Considering a dilution of 1:1 in PBS of bone marrow

		#Initialize the consumption of FP

		self.volume_FP_spent = 0

		#Insert the cost of Ficoll Paque per ml in EUR

		self.cost_FP = 0.39

		self.cost_pbs = 0.075

		#Insert the time in centrifuge in days (it's 10 minutes)

		self.time_centrifuge = 10/(60*24)

class ExpansionDB(object):

	'''Class containing all the base information for expansion of MSCs'''
	def __init__(self):
	
		self.exp_time = 5

class DownstreamDB(object):

	'''Class containing all relevant information for DSP'''

	def __init__(self):

		self.wash_vr_time = 4/24

		self.vr_yield = triangular(0.75,0.8,0.85)

		self.ff_yield = 1

		#For simplicity, assumed that, in each wash, the culture medium used for washing is equal to the total volume out
		#When cells are resuspended after harvesting.

		self.number_washes = 2

		self.wash_vr_time = self.wash_vr_time/self.number_washes

		#Cell concentration for storage per ml - for the release, choose 20 million cells

		self.cell_concentration = 12.5e6

		#time for finishing the vial filling and storage

		self.finish_time = 2/24

		#self.finish_time = 2

		self.cryovial_vol = 2

		self.cryovials_spent = 0

		#Assuming it's the 2ml vials from Nunc in the laboratory

		self.cost_cryovial = 1.27

		#Assuming it's the Gibco recovery medium

		self.cryomedium = 2.68

		self.cryomedium_ratio = 0.5

		self.cryomedium_spent = 0

class QualityControlDB(object):

	'''Class containing all relevant information for quality control testing'''

	def __init__(self):

		self.pass_release_ratio = 1

		#Cost of release testing per donor. Assuming the 10k generally considered.

		self.cost_release_testing = 10000

		#Assume that the quality controls take 7 days (Because of mycoplasma assays, potency, etc)

		self.qc_time = 7

class Database(FacilityDB,IsolationDB,ExpansionDB,DownstreamDB,QualityControlDB):

	'''Class that inherits from all the database instances when appropriate'''

	def __init__(self):

		FacilityDB.__init__(self)

		IsolationDB.__init__(self)

		ExpansionDB.__init__(self)

		DownstreamDB.__init__(self)

		QualityControlDB.__init__(self)