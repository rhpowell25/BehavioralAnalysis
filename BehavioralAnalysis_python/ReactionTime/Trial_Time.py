# -*- coding: utf-8 -*-

from Plot_Specs import Font_Specs
from TrialLength__RxnTime import TrialLength_RxnTime
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

class Trial_Time_Data():
    
    def __init__(self, xds_morn, xds_noon):

        #%% Display the funtion being used
        print('Reaction Time Box Plot:')
        
        # Date
        File_Name = xds_noon._lab_data__meta['raw_file_name']
        split_info = File_Name.split('_')
        self.Date = split_info[0]
        # Task
        self.Task = split_info[2]
        # Drug
        if 'Caff' in File_Name:
            self.Drug = 'Caffeine'
        if 'Lex' in File_Name:
            self.Drug = 'Escitalopram'
        if 'Cyp' in File_Name:
            self.Drug = 'Cypro'
        if 'Con' in File_Name:
            self.Drug = 'Control'
        
        # Font & plotting specifications
        self.font_specs = Font_Specs()
        
        #%% Get the reaction times & trial variables
        
        print('Morning')
        Rxn_Time_Vars_morn = TrialLength_RxnTime(xds_morn)
        self.target_dirs_morn = Rxn_Time_Vars_morn.target_dirs
        self.unique_tgts_morn = Rxn_Time_Vars_morn.target_centers
        self.rxn_time_morn = Rxn_Time_Vars_morn.rxn_time
        self.trial_length_morn = Rxn_Time_Vars_morn.trial_length
        
        print('Afternoon')
        Rxn_Time_Vars_noon = TrialLength_RxnTime(xds_noon)
        self.target_dirs_noon = Rxn_Time_Vars_noon.target_dirs
        self.unique_tgts_noon = Rxn_Time_Vars_noon.target_centers
        self.rxn_time_noon = Rxn_Time_Vars_noon.rxn_time
        self.trial_length_noon = Rxn_Time_Vars_noon.trial_length
        
        #%% Check to see if both sessions use a consistent number of target centers
        
        if np.array_equal(self.unique_tgts_morn, self.unique_tgts_noon) == False:
            print('Uneven target centers between morning & afternoon')
            # Only use the info of target centers conserved between morn & noon
            shared_target_centers_idx = np.in1d(self.unique_tgts_morn, self.unique_tgts_noon)
            self.shared_target_centers = self.unique_tgts_morn[np.where(shared_target_centers_idx)[0][0]]
        else:
            self.shared_target_centers = self.unique_tgts_morn
            
         #%% Check to see if both sessions use a consistent number of target directions   
            
        if np.array_equal(self.target_dirs_morn, self.target_dirs_noon) == False:
            print('Uneven target directions between morning & afternoon')
            return
        else:
            self.shared_target_dirs = self.target_dirs_morn
            
     
    #%% Reaction Time Box Plots
    
    def Rxn_Time_BoxPlot(self, Save_Figs):
    
        #%% Put the trial times in the same struct of the plots

        rxn_time = [[] for ii in range(len(self.rxn_time_morn))]
        for pp in range(len(self.rxn_time_morn)):
            rxn_time[pp].append(self.rxn_time_morn[pp])
            rxn_time[pp].append(self.rxn_time_noon[pp])
            
        #%% Do the statistics
        
        rxn_time_p_val = np.zeros(len(self.rxn_time_morn))
        for pp in range(len(self.rxn_time_morn)):
            rxn_time_p_val[pp] = stats.ttest_ind(self.rxn_time_morn[pp], self.rxn_time_noon[pp])[1]
            
        #%% Find the percent change
        
        avg_rxn_time_morn = np.zeros(len(self.rxn_time_morn))
        avg_rxn_time_noon = np.zeros(len(self.rxn_time_noon))
        rxn_time_perc_change = np.zeros(len(self.rxn_time_morn))
        for pp in range(len(self.rxn_time_morn)):
            avg_rxn_time_morn[pp] = np.mean(self.rxn_time_morn[pp])
            avg_rxn_time_noon[pp] = np.mean(self.rxn_time_noon[pp])
            rxn_time_perc_change[pp] = avg_rxn_time_noon[pp] - avg_rxn_time_morn[pp] / abs(avg_rxn_time_morn[pp])
            
        #%% Find the y-axis limits
        
        rxn_time_y_min = np.zeros(len(self.rxn_time_morn))
        rxn_time_y_max = np.zeros(len(self.rxn_time_morn))

        #%% Plot the reaction time data
    
        for pp in range(len(self.rxn_time_morn)):
            
            # Put in the categories for the x-axis
            x_axis = np.array([1, 2])
            categ = ['Morning', 'Afternoon']
            
            # Boxplot
            rxn_time_fig, figure_axes = plt.subplots()
            plt.boxplot(rxn_time[pp])
            plt.xticks(x_axis, categ)
            
            # Increase the axes font
            plt.xticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
            plt.yticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5) 
            
            # Collect the current y-axis limits
            y_limits = figure_axes.get_ylim()
            rxn_time_y_min[pp] = y_limits[0]
            rxn_time_y_max[pp] = y_limits[1]
            
            # Annotation of the n-count
            morn_succ_trials = 'n = ' + str(len(rxn_time[pp][0]))
            noon_succ_trials = 'n = ' + str(len(rxn_time[pp][1]))
            
            plt.text(0.35, 0.05, morn_succ_trials, verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
            plt.text(0.65, 0.05, noon_succ_trials, verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
    
            # Annotation of the p-value
            if round(rxn_time_p_val[pp], 3) > 0:
                plt.text(0.85, 0.9, 'p = ' + str(round(rxn_time_p_val[pp], 3)), 
                         verticalalignment = 'center', horizontalalignment = 'center', 
                         transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
            if round(rxn_time_p_val[pp], 3) == 0:
                plt.text(0.85, 0.9, 'p < 0.001', verticalalignment = 'center', horizontalalignment = 'center', 
                         transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
     
            # Annotation of the percent change
            if round(rxn_time_perc_change[pp], 3) != 0:
                plt.text(0.15, 0.9, 'PC = ' + str(round(rxn_time_perc_change[pp], 3)), 
                         verticalalignment = 'center', horizontalalignment = 'center', 
                         transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
            if round(rxn_time_perc_change[pp], 3) == 0:
                plt.text(0.15, 0.9, 'PC ≈ 0', verticalalignment = 'center', horizontalalignment = 'center', 
                         transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
    
            # Title
            if self.Task != 'PG':
                title_string = self.Date + ' ' + self.Task + ' ' + self.Drug + ': ' + str(int(self.shared_target_dirs[pp])) + '° ' + 'TgtCenter at ' + str(self.shared_target_centers[pp])
            else:
                title_string = self.Date + ' ' + self.Task + ' ' + self.Drug + ': ' + 'TgtCenter at ' + str(self.shared_target_centers[pp])
    
            plt.title(title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
    
            # Labels
            plt.ylabel('Reaction Time (Sec.)', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
            # Axis limits
            plt.xlim([0, 3])
            plt.ylim(rxn_time_y_min[pp] - 0.06, rxn_time_y_max[pp])
            
            # Figure Saving
            if Save_Figs != 0:
                save_dir = 'C:/Users/rhpow/Desktop/'
                fig_title = title_string
                fig_title = str.replace(fig_title, ':', '')
                fig_title = str.replace(fig_title, 'vs.', 'vs')
                fig_title = str.replace(fig_title, 'mg.', 'mg')
                fig_title = str.replace(fig_title, 'kg.', 'kg')
                fig_title = str.replace(fig_title, '.', '_')
                fig_title = str.replace(fig_title, '/', '_')
                plt.savefig(save_dir + fig_title + '.' + Save_Figs)
                plt.close()
        
    #%% Trial Length Box Plots
    
    def Trial_Length_BoxPlot(self, Save_Figs):
        
        #%% Put the trial times in the same struct of the plots

        trial_length = [[] for ii in range(len(self.trial_length_morn))]
        for pp in range(len(self.trial_length_morn)):
            trial_length[pp].append(self.trial_length_morn[pp])
            trial_length[pp].append(self.trial_length_noon[pp])
        
        #%% Do the statistics
     
        trial_length_p_val = np.zeros(len(self.trial_length_morn))
        for pp in range(len(self.trial_length_morn)):
            trial_length_p_val[pp] = stats.ttest_ind(self.trial_length_morn[pp], self.trial_length_noon[pp])[1]
    
        #%% Find the percent change
    
        avg_trial_length_morn = np.zeros(len(self.trial_length_morn))
        avg_trial_length_noon = np.zeros(len(self.trial_length_noon))
        trial_length_perc_change = np.zeros(len(self.trial_length_morn))
        for pp in range(len(self.trial_length_morn)):
            avg_trial_length_morn[pp] = np.mean(self.trial_length_morn[pp])
            avg_trial_length_noon[pp] = np.mean(self.trial_length_noon[pp])
            trial_length_perc_change[pp] = avg_trial_length_noon[pp] - avg_trial_length_morn[pp] / abs(avg_trial_length_morn[pp])
    
        #%% Find the y-axis limits

        trial_length_y_min = np.zeros(len(self.trial_length_morn))
        trial_length_y_max = np.zeros(len(self.trial_length_morn))
    
        #%% Plot the trial length data
    
        for pp in range(len(self.trial_length_morn)):
            
            # Put in the categories for the x-axis
            x_axis = np.array([1, 2])
            categ = ['Morning', 'Afternoon']
            
            # Boxplot
            rxn_time_fig, figure_axes = plt.subplots()
            plt.boxplot(trial_length[pp])
            plt.xticks(x_axis, categ)
            
            # Increase the axes font
            plt.xticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
            plt.yticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5) 
            
            # Collect the current y-axis limits
            y_limits = figure_axes.get_ylim()
            trial_length_y_min[pp] = y_limits[0]
            trial_length_y_max[pp] = y_limits[1]
            
            # Annotation of the n-count
            morn_succ_trials = 'n = ' + str(len(trial_length[pp][0]))
            noon_succ_trials = 'n = ' + str(len(trial_length[pp][1]))
            
            plt.text(0.35, 0.05, morn_succ_trials, verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
            plt.text(0.65, 0.05, noon_succ_trials, verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
    
            # Annotation of the p-value
            if round(trial_length_p_val[pp], 3) > 0:
                plt.text(0.85, 0.9, 'p = ' + str(round(trial_length_p_val[pp], 3)), 
                         verticalalignment = 'center', horizontalalignment = 'center', 
                         transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
            if round(trial_length_p_val[pp], 3) == 0:
                plt.text(0.85, 0.9, 'p < 0.001', verticalalignment = 'center', horizontalalignment = 'center', 
                         transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
     
            # Annotation of the percent change
            if round(trial_length_perc_change[pp], 3) != 0:
                plt.text(0.15, 0.9, 'PC = ' + str(round(trial_length_perc_change[pp], 3)), 
                         verticalalignment = 'center', horizontalalignment = 'center', 
                         transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
            if round(trial_length_perc_change[pp], 3) == 0:
                plt.text(0.15, 0.9, 'PC ≈ 0', verticalalignment = 'center', horizontalalignment = 'center', 
                         transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
    
            # Title
            if self.Task != 'PG':
                title_string = self.Date + ' ' + self.Task + ' ' + self.Drug + ': ' + str(int(self.shared_target_dirs[pp])) + '° ' + 'TgtCenter at ' + str(self.shared_target_centers[pp])
            else:
                title_string = self.Date + ' ' + self.Task + ' ' + self.Drug + ': ' + 'TgtCenter at ' + str(self.shared_target_centers[pp])
    
            plt.title(title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
    
            # Labels
            plt.ylabel('Trial Length (Sec.)', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
            # Axis limits
            plt.xlim([0, 3])
            plt.ylim(trial_length_y_min[pp] - 0.06, trial_length_y_max[pp])
            
            # Figure Saving
            if Save_Figs != 0:
                save_dir = 'C:/Users/rhpow/Desktop/'
                fig_title = title_string
                fig_title = str.replace(fig_title, ':', '')
                fig_title = str.replace(fig_title, 'vs.', 'vs')
                fig_title = str.replace(fig_title, 'mg.', 'mg')
                fig_title = str.replace(fig_title, 'kg.', 'kg')
                fig_title = str.replace(fig_title, '.', '_')
                fig_title = str.replace(fig_title, '/', '_')
                plt.savefig(save_dir + fig_title + '.' + Save_Figs)
                plt.close()



    
    