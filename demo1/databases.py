# databases.py - contain all the classes associated with fixed values

import math
import pandas as pd
import numpy as np
from numpy.random import *

# FacilityDB - Class containing variables to pass to stemfactory.py


class FacilityDB(object):

    def __init__(self):

        # FACILITY STRUCTURE PARAMETERS

        self.facility_area = 400  # in sq mt

        self.clean_room_ratio = 0.2

        self.cost_sqmt_clean_room = 5815

        self.cost_sqmt_other_room = 3392

        total_cost_clean_room = (self.cost_sqmt_clean_room
                                    * self.clean_room_ratio)

        total_cost_other_room = (self.cost_sqmt_other_room
                                    * (1 - self.clean_room_ratio))

        self.total_facility_cost = (self.facility_area
                                    * (total_cost_clean_room
                                        + total_cost_other_room))

        self.facility_depreciation_period = 15*365.25  # Period in days

        co2_daily_cost = 6000/365.25

        gases_daily_cost = 15600/365.25

        add_supplies_daily_cost = 7900/365.25

        requal_daily_cost = 65400/365.25

        maint_daily_cost = 52800/365.25

        cleaning_daily_cost = 28000/365.25

        garment_daily_cost = 2000/365.25

        compilation_daily_cost = (co2_daily_cost + gases_daily_cost
                                + add_supplies_daily_cost + requal_daily_cost
                                + maint_daily_cost + cleaning_daily_cost
                                + garment_daily_cost)

        self.daily_facility_cost = (self.total_facility_cost
                                / self.facility_depreciation_period
                                + compilation_daily_cost)

        self.total_workers = 3

        self.total_bscs = 4

        self.total_incubators = 8

        self.total_centrifuges = 4

        self.occupied_inc = 0  # Provides the initial no. occupied incubators

        self.daily_worker = 250  # Daily rate of each worker

        self.unit_incubator_cost = 16000

        self.unit_bsc_cost = 15300

        self.unit_centrifuge_cost = 12000

        self.equipment_depreciation_period = 5*365.25  # In days

        self.daily_incubator_cost = (self.unit_incubator_cost
                                    * self.total_incubators
                                    / self.equipment_depreciation_period)

        self.daily_bsc_cost = (self.unit_bsc_cost
                                * self.total_bscs
                                / self.equipment_depreciation_period)

        self.daily_centrifuge_cost = (self.unit_centrifuge_cost
                                    * self.total_centrifuges
                                    / self.equipment_depreciation_period)

        # DATABASE OF PLANAR TECHNOLOGIES

        self.names_of_planar_ets = ['t-flask25',
                       't-flask75',
                       't-flask175',
                       't-flask225',
                       'cellstack1',
                       'cellstack2',
                       'cellstack5',
                       'cellstack10']

        self.area_planar_ets = pd.Series(data=([25, 75, 175, 225, 636,
                                636*2, 636*5, 636*10]),
                                index=self.names_of_planar_ets)  # In sq. cm

        self.media_volume_planar_ets = pd.Series(data=([5, 15, 35, 45,
                                    130, 260, 650, 1300]),
                                    index=self.names_of_planar_ets)  # In mL

        self.harvesting_volume_planar_ets = pd.Series(data=([1.75, 3.5, 7,
                                            9, 25, 50, 125, 250]),
                                            index=self.names_of_planar_ets)

        self.units_resource_planar = pd.DataFrame(
                                    data=[[100, 1],
                                            [100, 1],
                                            [100, 1],
                                            [100, 1],
                                            [60, 1],
                                            [60, 1],
                                            [24, 1],
                                            [12, 1]],
                                    index=self.names_of_planar_ets,
                                    columns=['incubator', 'worker']
                                    )  # No. of planar flasks

        self.unit_planar_ets_cost = pd.Series(data=[1.65, 4.85, 7.38, 8.55,
                                            33.89, 60.12, 131.96, 142.72],
                                            index=self.names_of_planar_ets)

        self.et_operation_times_planar = pd.DataFrame(
            data=[[0.38/(24*10), 0.38/(24*10), 0.5/(24*10)],
                   [0.38/(24*10), 0.38/(24*10), 0.5/(24*10)],
                   [0.38/(24*10), 0.38/(24*10), 0.5/(24*10)],
                   [0.38/(24*10), 0.38/(24*10), 0.5/(24*10)],
                   [0.15/24, 0.15/24, 0.06/24],
                   [0.15/24, 0.15/24, 0.06/24],
                   [0.20/24, 0.20/24, 0.16/24],
                   [0.25/24, 0.25/24, 0.25/24]],
            index=self.names_of_planar_ets,
            columns=['seeding', 'feeding', 'harvesting']
            )

        self.harvesting_time_incubator = 14/(60*24)

        self.time_bw_feedings = 3

        # COSTS OF COMMON REAGENTS

        self.cost_dmem_fbs = 0.04*0.90 + 2.12*0.10  # FBS

        self.cost_basal_medium = 0.04

        self.cost_tryple = 0.21

        # Initialize the volumes spent for each reagent

        self.volume_media_spent = 0

        self.volume_basal_medium_spent = 0

        self.volume_tryple_spent = 0

        self.planar_ets_spent = pd.Series(data=([0]
                                * len(self.names_of_planar_ets)),
                                index=self.names_of_planar_ets)

        self.total_doses = 0

        self.processed_donors = 0


class IsolationDB(object):

    # IsolationDB - Isolation unit operation specific parameters

    def __init__(self):

        self.volume_BM = 10  # in mL

        self.volume_FP = self.volume_BM * 2  # 1:1 dilution in PBS

        self.volume_FP_spent = 0  # Initialize BM consumption

        self.cost_FP = 0.39  # Cost/mL

        self.cost_pbs = 0.075  # Cost/mL

        self.time_centrifuge = 10/(60*24)  # Time centrifuge in days


class ExpansionDB(object):

    # ExpansionDB - Expansion unit operation specific parameters

    def __init__(self):

        self.exp_time = 5

        self.max_passages_exp = 5 

        self.confluence_passage_exp = [7]*(self.max_passages_exp+1)

        self.total_expansion_time = sum(self.confluence_passage_exp)


class DifferentiationDB(object):

    # DifferentiationDB - Differentiation unit operation specific parameters

    def __init__(self):

        self.no_diff_stages = 6  # Stages = no. intermediate cell types

        self.diff_stages = ['S1',
                       'S2',
                       'S3',
                       'S4',
                       'S5',
                       'S6']

        # Differentiation requires an initial cell aggregation

        self.cost_aggregation_media_ml = 0  # From StemPro XF ESCs medium

        self.time_aggregation = 0  # in days

        self.cost_differentiation_media_ml = pd.Series(data=[4.32, 2.31, 2.40,
                                                2.38, 1.16, 0.21],
                                                index=self.diff_stages)

        # self.time_stage annotates the time points of completion
        # in each differentiation stage

        self.time_stage = pd.Series(data=[3, 6, 8, 13, 20, 30],
                                    index=self.diff_stages)

        # self.number_feedings_stage annotates how many culture medium
        # exchanges (feedings) are performed per differentiation stage
        # to be computer by the donor.py and expansiontechnology.py
        # modules

        self.number_feedings_stage = pd.Series(data=[2, 2, 2, 3, 7, 5],
                                                index=self.diff_stages)

        self.total_differentiation_time = self.time_stage[-1]


class DownstreamDB(object):

    # DownstreamDB - DSP unit operation specific parameters

    def __init__(self):

        self.wash_vr_time = 4/24

        self.vr_yield = 0.8  # Ratio cells out after volume reduction

        self.ff_yield = 1  # Ratio cells out after fill finish

        self.number_washes = 2  # Total no. centrifugation cycles

        self.wash_vr_time = self.wash_vr_time/self.number_washes

        self.cell_concentration = 12.5e6  # Cells/mL for storage

        self.finish_time = 2/24  # time for finishing the vial filling and storage

        self.cryovial_vol = 2  # in mL

        self.cryovials_spent = 0  # initialize the storage

        self.cost_cryovial = 1.27  # unit price

        self.cryomedium = 2.68  # cost/mL of Gibco recovery medium

        self.cryomedium_ratio = 0.5  # ratio cryomedium to basal medium

        self.cryomedium_spent = 0  # initialize volume cryomedium

        self.total_dsp_time = self.wash_vr_time + self.finish_time


class QualityControlDB(object):

    # QualityControlDB - QC unit operation specific parameters

    def __init__(self):

        self.pass_release_ratio = 0.9  # Ratio final batches passing QC

        # Cost of release testing per donor.
        # Assuming the 10k generally considered for all types of testing
        # Identity, Potency, Sterility, Purity

        self.cost_release_testing = 10000

        self.qc_time = 3/24  # 3h of QC for fresh infusion

class Database(FacilityDB, IsolationDB, ExpansionDB,
                DifferentiationDB, DownstreamDB, QualityControlDB):

    # Incorporate all the other classes

    def __init__(self):

        FacilityDB.__init__(self)

        IsolationDB.__init__(self)

        ExpansionDB.__init__(self)

        DifferentiationDB.__init__(self)

        DownstreamDB.__init__(self)

        QualityControlDB.__init__(self)


        # Incorporates the annual demands

        self.annual_cell_demand = 10e9

        self.cells_per_batch = 1e9
    
        self.cells_per_dose = 10e6

        self.doses_per_batch = round(self.cells_per_batch / self.cells_per_dose,0)

        self.annual_batches = math.ceil((self.annual_cell_demand / self.cells_per_batch)
                                / self.pass_release_ratio)

        # Updates the dimensioning of the campaign as such

        self.total_campaign_time = self.total_expansion_time + self.total_dsp_time + self.qc_time

        self.annual_campaigns = math.ceil(365.25 / self.total_campaign_time)

        # Since this is a scale out process, assumes that each incubator will fill the needs out

        self.total_incubators = math.ceil(self.annual_batches / self.annual_campaigns)

        self.total_bscs = self.total_incubators / 2

        self.total_centrifuges = self.total_incubators / 2

        self.total_workers = 3 * round(self.total_incubators / 2, 0)

        # Updates the daily costs of equipments

        self.daily_incubator_cost = (self.unit_incubator_cost
                                    * self.total_incubators
                                    / self.equipment_depreciation_period)

        self.daily_bsc_cost = (self.unit_bsc_cost
                                * self.total_bscs
                                / self.equipment_depreciation_period)

        self.daily_centrifuge_cost = (self.unit_centrifuge_cost
                                    * self.total_centrifuges
                                    / self.equipment_depreciation_period)

        # Needs a debate on the best method, installed capacity or fixed capacity

