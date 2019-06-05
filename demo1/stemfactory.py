# stemfactory.py - Class responsible for keeping the facility working

import os
import sys
import simpy
import random
import math
import logging
import time
import pandas as pd

import matplotlib.pyplot as plt
import pylab as p
import seaborn as sns

# Call the other custom modules

from donor import *

from databases import *


class StemFactory():

    def __init__(self, env, db):

        self.env = env  # Necessary parameter for discrete event simulation

        # INITIALIZING COSTS

        self.cost_ET = 0  # Consumable costs (ET = Expansion Technologies)

        self.cost_reagents = 0

        self.cost_building = 0

        self.cost_equipment = 0

        self.cost_labor = 0

        self.cost_qc = 0

        self.total_costs = 0

        self.costs_per_dose = 0

        self.total_doses = 0  # Initialize the total doses produced

        self.number_donors = db.annual_batches  # Maximum number of donors/batches to process

        self.processed_donors = 0  # Initialize count of processed donors

        # INITIALIZING QUEUES FOR EQUIPMENT

        self.inc_list = simpy.Resource(env, db.total_incubators)

        self.bsc_list = simpy.Resource(env, db.total_bscs)

        # INITIALIZING FRAMES FOR STORAGE OF COSTS AND TIMES

        self.cost_donor_stage = pd.DataFrame(
                                    data=np.zeros((self.number_donors,5)),
                                    index=range(self.number_donors),
                                    columns=['isolation','expansion','differentiation','dsp','release']
                                )

        self.cost_donor_variable = pd.DataFrame(
                                    data=np.zeros((self.number_donors, 6)),
                                    index=range(self.number_donors),
                                    columns=['consumables', 'reagents',
                                            'qc', 'building',
                                            'equipment', 'labor']
                                    )

        self.costs_resources_plot = pd.DataFrame(data=np.zeros((self.number_donors*6,2)),
                                    index=range(self.number_donors*6),
                                    columns=['resources','cost'])

        self.times_process_stage = pd.DataFrame(
                                    data = np.zeros((self.number_donors,5)),
                                    index = range(self.number_donors),
                                    columns = ['isolation','expansion','differentiation','dsp','release']
                                    )

        self.costs_stages_plot = pd.DataFrame(data=np.zeros((self.number_donors*5,2)),
                                    index=range(self.number_donors*5),
                                    columns=['stages','cost'])

        self.time_processing_histogram = pd.Series(
            data = np.zeros((self.number_donors)),
            index = range(self.number_donors)
            )

        self.total_costs_dose = pd.Series(
            data=np.zeros((self.number_donors)),
            index=range(self.number_donors)
            )

        self.total_doses_per_donor = pd.Series(
            data=np.zeros((self.number_donors)),
            index=range(self.number_donors)
            )

        self.total_donor_costs_histogram = pd.Series(
            data=np.zeros((self.number_donors)),
            index=range(self.number_donors)
            )


        self.costs_failed_doses = 0  # Initialize for redistribution


    def run(self, env, db):

        # run is the method coordinating how the donors are processed

        for donor_index in range(self.number_donors):

            # Launches each donor

            env.process(self.donor(env, donor_index, db))

        while db.processed_donors < self.number_donors:

            # Time keeps on advancing while donors are not complete

            yield env.timeout(0.1)


        #print('Donor costs in histogram before fixed cost calulator %.2f' % self.total_donor_costs_histogram)

    
        # Calculate the costs

        self.fixed_cost_calculator(env, db)

        #print('Donor costs in histogram after fixed cost calulator %.2f' % self.total_donor_costs_histogram)
        # Later, insert the method for storing the CSV files of each relevant histogram

        self.csv_generator(env, db)

        # Finally, generates the plots automatically

        self.plot_generator(env, db)

    def donor(self, env, donor_index, db):

        # Create a donor from the class

        dn = Donor(env, donor_index, db)

        with self.inc_list.request() as request:

            # Donors are waiting in the queue for incubator space

            yield request

            # As a donor enters the incubator, its culture begins

            donor_proc = dn.processing(env, db)

            yield env.process(donor_proc)

            while dn.processed == 0:

                # Time advances while donor is not fully processed

                yield env.timeout(0.0001)

            # Sum the number of doses produced

            self.total_doses += dn.doses_no

            # Update the processing time

            self.time_processing_histogram[donor_index] = dn.time_processing
            
            # Adds the variable costs of the donor to the dataframes

            self.variable_cost_calculator(env, dn, db)

            # Add doses to the number of final doses if they passed the quality controls

            if dn.doses_no > 0:

                self.total_doses_per_donor[dn.donor_index] = dn.doses_no

            # Add the consumables spent per donor
            # Correct this later

            
            # for et in db.names_of_planar_ets:

               # self.ets_per_donor.loc[dn.donor_index,et] = dn.planar_ets_spent.loc[et]


    def variable_cost_calculator(self, env, dn, db):

        # Method for compilation of all costs that are donor dependent

        # Add the costs per stage

        self.cost_donor_stage.loc[dn.donor_index, 'expansion'] = dn.total_costs_expansion

        self.cost_donor_stage.loc[dn.donor_index, 'differentiation'] = dn.total_costs_differentiation

        self.cost_donor_stage.loc[dn.donor_index, 'dsp'] = dn.total_costs_DSP

        self.cost_donor_stage.loc[dn.donor_index, 'release'] = dn.total_costs_release

        # Add the costs per variable

        self.cost_donor_variable.loc[dn.donor_index,'consumables'] = dn.total_costs_ets

        self.cost_donor_variable.loc[dn.donor_index,'reagents'] = dn.total_costs_reagents

        self.cost_donor_variable.loc[dn.donor_index,'qc'] = dn.total_costs_release

        # Add the total variable costs to the histogram

        self.total_donor_costs_histogram[dn.donor_index] = dn.total_costs_donor

        #print('Total costs of donor %.2f' % dn.total_costs_donor)

        #print('Total costs of donor histogram %.2f' % self.total_donor_costs_histogram[dn.donor_index])

    def fixed_cost_calculator(self, env, db):

        # Add all the fixed costs that vary with the total number of days in culture

        total_process_time = sum(self.time_processing_histogram)

        #print('Total process time inside fixed cost calculator %d' % total_process_time)

        self.cost_building = db.daily_facility_cost * math.ceil(env.now)

        self.cost_equipment = ((db.daily_incubator_cost * db.total_incubators
                                + db.daily_bsc_cost * db.total_bscs
                                + db.daily_centrifuge_cost * db.total_centrifuges)
                                * math.ceil(env.now))

        self.cost_labor = ((db.daily_worker * db.total_workers)
                            * math.ceil(env.now))

        #Sum these costs to the total costs of each donor

        for i in range(self.number_donors):

            self.cost_building_per_donor = (self.cost_building 
                                            * (self.time_processing_histogram[i]
                                                / total_process_time))

            self.cost_equipment_per_donor = (self.cost_equipment
                                            * (self.time_processing_histogram[i]
                                                / total_process_time))

            self.cost_labor_per_donor = (self.cost_labor
                                        * (self.time_processing_histogram[i]
                                            / total_process_time))

            self.cost_donor_variable.loc[i,'building'] = self.cost_building_per_donor

            self.cost_donor_variable.loc[i,'equipment'] = self.cost_equipment_per_donor

            self.cost_donor_variable.loc[i,'labor'] = self.cost_labor_per_donor

            self.total_donor_costs_histogram[i] += (self.cost_building_per_donor
                                                    + self.cost_equipment_per_donor
                                                    + self.cost_labor_per_donor)

            if self.total_doses_per_donor[i] == 0:

                self.costs_failed_doses += self.total_donor_costs_histogram[i]

            else:

                # Divide the costs per donor to get the doses

                self.total_costs_dose[i] = (self.total_donor_costs_histogram[i]
                    / self.total_doses_per_donor[i])

            # Add the fixed costs
            # FIX THIS SO THAT IT CAN GET WITH BETTER LINE NUMBERS

            self.cost_donor_stage.loc[i,'expansion'] += (self.cost_building_per_donor + self.cost_equipment_per_donor+self.cost_labor_per_donor)*round(self.times_process_stage.loc[i,'expansion']/self.time_processing_histogram[i],2)

            self.cost_donor_stage.loc[i,'differentiation'] += (self.cost_building_per_donor+self.cost_equipment_per_donor+self.cost_labor_per_donor)*round(self.times_process_stage.loc[i,'differentiation']/self.time_processing_histogram[i],2)
            self.cost_donor_stage.loc[i,'dsp'] += (self.cost_building_per_donor+self.cost_equipment_per_donor+self.cost_labor_per_donor)*round(self.times_process_stage.loc[i,'dsp']/self.time_processing_histogram[i],2)

            self.cost_donor_stage.loc[i,'release'] += (self.cost_building_per_donor+self.cost_equipment_per_donor+self.cost_labor_per_donor)*round(self.times_process_stage.loc[i,'release']/self.time_processing_histogram[i],2)

        #Add the costs of the failed doses equally distributed for all

        for i in range(self.number_donors):

            if self.total_doses_per_donor[i] > 0:

                self.total_costs_dose[i] += self.costs_failed_doses/sum(self.total_doses_per_donor)

        #Ended

    def csv_generator(self, env, db):
        
        # Create the CSV commands for each of the platforms

        # Creates a new folder to store the files

        cwd = os.getcwd()

        new_dir = cwd + '/csv_files'

        os.makedirs(new_dir)

        os.chdir(new_dir)

        self.time_processing_histogram.to_csv('process_times.csv')

        self.cost_donor_variable.to_csv('cost_donor_variable.csv')

        self.cost_donor_stage.to_csv('cost_donor_stage.csv')

        #print('Donor costs in histogram csv generator %.2f' % self.total_donor_costs_histogram)

        self.total_donor_costs_histogram.to_csv('total_costs_donor.csv')

        self.total_doses_per_donor.to_csv('total_doses_donor.csv')

        self.total_costs_dose.to_csv('total_costs_dose.csv')

        # Changes back to the root directory

        os.chdir(cwd) 


    def plot_generator(self, env, db):

        # Create the commands for plots to be stored

        # Start by iterating over the total number of donors

        for dn_id in range(self.number_donors):

            # Resources

            self.costs_resources_plot.loc[dn_id*6:dn_id*6+5,'resources'] = ['consumables', 'reagents',
                                                                'qc', 'building',
                                                                'equipment', 'labor']

            self.costs_resources_plot.loc[dn_id*6:dn_id*6+5,'cost'] = np.array(self.cost_donor_variable.iloc[dn_id,:]) 
        
            # Stages

            self.costs_stages_plot.loc[dn_id*5:dn_id*5+4,'stages'] = ['isolation','expansion','differentiation','dsp','release']

            self.costs_stages_plot.loc[dn_id*5:dn_id*5+4,'cost'] = np.array(self.cost_donor_stage.iloc[dn_id,:]) 

        # Create a new directory to store the plots

        cwd = os.getcwd()

        new_dir = cwd + '/plots'

        os.makedirs(new_dir)

        os.chdir(new_dir)

        # Do the plots

        self.barplot_parameters('cost', 'resources', self.costs_resources_plot, 'CostsResources')

        self.barplot_parameters('cost', 'stages', self.costs_stages_plot, 'CostsStages')

        if self.number_donors > 1:

            # Only makes sense to generate the histogram if there is a variability

            self.histplot_parameters('cost', 'dose', self.total_costs_dose.to_frame(), 'CostsPerDose')

            self.histplot_parameters('doses', 'donor', self.total_doses_per_donor.to_frame(), 'DosesPerDonor')

        # Change back to the main directory

        os.chdir(cwd)


    def barplot_parameters(self, xlabel, ylabel, data_frame, plot_name):

        # Defines the Seaborn style for research manuscript

        sns.set(color_codes=True)
        sns.set_style("white")
        sns.set_style("ticks")
        sns.set_context("paper")

        pal = sns.color_palette()
        pal.as_hex()

        # Create the Matplotlib/Seaborn plot handle

        ax = sns.barplot(x=xlabel, y=ylabel, data=data_frame, color="b")

        fig = ax.get_figure()

        # Define the size for a one column panel

        width = 3.2
        height = 2.4

        # Make an appropriate scale

        ax.autoscale(tight=False)
        ax.margins(0.05)
        ax.grid(False)

        string_for_title = xlabel.title() + ' per ' + ylabel.title()

        plt.title(string_for_title)
        plt.xlabel('Cost ($)')
        plt.ylabel(' ')
        fig.set_size_inches(width, height)

        string_for_name = plot_name + '.tif'

        fig.savefig(string_for_name,dpi=600,bbox_inches='tight')

        plt.close()

    def histplot_parameters(self, xlabel, ylabel, data_frame, plot_name):

        # Defines the Seaborn style for research manuscript

        sns.set(color_codes=True)
        sns.set_style("white")
        sns.set_style("ticks")
        sns.set_context("paper")

        pal = sns.color_palette()
        pal.as_hex()

        # Create the Matplotlib/Seaborn plot handle

        ax = sns.distplot(data_frame, kde=False, bins=100)

        fig = ax.get_figure()

        # Define the size for a one column panel

        width = 3.2
        height = 2.4

        # Make an appropriate scale

        ax.autoscale(tight=False)
        ax.margins(0.05)
        ax.grid(False)

        string_for_title = xlabel.title() + ' per ' + ylabel.title()

        plt.title(string_for_title)
        plt.xlabel(xlabel)
        plt.ylabel(' ')
        fig.set_size_inches(width, height)

        string_for_name = plot_name + '.tif'

        fig.savefig(string_for_name,dpi=600,bbox_inches='tight')

        plt.close()
