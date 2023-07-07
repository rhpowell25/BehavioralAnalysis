#%% Import basic packages
import numpy as np
import matplotlib.pyplot as plt
from Plot_Specs import Font_Specs
import matplotlib.font_manager as fm
    
#%% File Description:

# Plots the average cursor position, velocity, or acceleration of two 
# XDS files overtop one another
#
# -- Inputs --
# xds_morn: the first xds file
# xds_noon: the second xds file
# signal_choice: 'Pos', 'Vel', or 'Acc'
# event: 'trial_gocue', 'window_trial_gocue', 'trial_end', 
# 'window_trial_end', 'force_onset', 'window_force_onset', 'force_max', 
# 'window_force_max', 'window_force_deriv', 'force_deriv', 'cursor_onset', 
# 'window_cursor_onset', 'cursor_veloc', 'window_cursor_veloc', 
# 'cursor_acc', 'window_cursor_acc', 'EMG_max', 'window_EMG_max', 
# 'task_onset', or 'window_task_onset'
# Save_Figs: 'pdf', 'png', 'fig', or 0

def OverlapPlotCursor(xds_morn, xds_noon, signal_choice, event, Save_Figs):
    
    #%% Display the function being used
    print('Overlap Cursor Function:')

    #%% Collect the Y-Limits
    from CursorYLimit import CursorYLimit
    cursor_YLims = CursorYLimit(xds_morn, xds_noon, signal_choice)
    
    # End the function if there is no Y-Limit
    if cursor_YLims == 'NaN':
        print("There is no Y-Limit")
        return

    #%% Basic settings, some variable extractions, & definitions
    
    # Event lengths
    before_event = 3
    after_event = 3
                    
    if 'gocue' in event or 'force_onset' in event:
        # Define the window for the baseline phase
        time_before_gocue = 0.4
    elif 'end' in event:
        # Define the window for the movement phase
        time_before_end = xds_morn._lab_data__meta['TgtHold']
        
    # Font & plotting specifications
    font_specs = Font_Specs()
                    
    #%% Load the average cursor signals
    from AvgCursorSignal import AvgCursorSignal
    avg_curs_sig_morn = AvgCursorSignal(xds_morn, signal_choice, event, 'NaN', cursor_YLims, 0, 0)
    avg_curs_sig_noon = AvgCursorSignal(xds_noon, signal_choice, event, 'NaN', cursor_YLims, 0, 0)
                    
    #%% Extract the target directions & centers
    from Identify_Targets import Identify_Targets
    target_vars_morn = Identify_Targets(xds_morn)
    target_dirs_morn = target_vars_morn.target_dirs
    target_centers_morn = target_vars_morn.target_centers
    target_vars_noon = Identify_Targets(xds_noon)
    target_dirs_noon = target_vars_noon.target_dirs
    target_centers_noon = target_vars_noon.target_centers
    
    #%% Check to see if both sessions use a consistent number of targets
    
    # Find matching targets between the two sessions
    from Match_Targets import Match_Targets
    Matching_Idxs = Match_Targets(target_dirs_morn, target_dirs_noon, \
                                      target_centers_morn, target_centers_noon)
    Matching_Idxs_Morn = Matching_Idxs.Matching_Idxs_Morn
    Matching_Idxs_Noon = Matching_Idxs.Matching_Idxs_Noon
    
    # Only use the info of target centers conserved between morn & noon
    if not all(ii is True for ii in Matching_Idxs_Morn) or not all(ii is True for ii in Matching_Idxs_Noon):
        print('Uneven Targets Between Morning & Afternoon')
        target_centers_morn = target_centers_morn[Matching_Idxs_Morn]
        target_dirs_morn = target_dirs_morn[Matching_Idxs_Morn]
        target_centers_noon = target_centers_noon[Matching_Idxs_Noon]
        target_dirs_noon = target_dirs_noon[Matching_Idxs_Noon]
        
    #%% X-axis
    cursor_time = np.linspace(-before_event, after_event, len(avg_curs_sig_morn[0]))
        
    #%% Y-axis limits 
    y_limits = np.zeros(2)
    y_max = np.zeros((len(avg_curs_sig_morn), 2))
    y_min = np.zeros((len(avg_curs_sig_morn), 2))
    
    for ii in range(len(avg_curs_sig_morn)):
        y_max[ii,0] = np.max(avg_curs_sig_morn[ii])
        y_max[ii,1] = np.max(avg_curs_sig_noon[ii])
        y_min[ii,0] = np.min(avg_curs_sig_morn[ii])
        y_min[ii,1] = np.min(avg_curs_sig_noon[ii])
        
    y_limits[0] = np.min(y_min) - 0.125
    y_limits[1] = np.max(y_max) + 0.125
        
    #%% Plot the two overlapped
    for ii in range(len(avg_curs_sig_morn)):
        
        
        fig, ax = plt.subplots(figsize = (11, 4))
        
        plt.plot(cursor_time, avg_curs_sig_morn[ii], color = (0.9290, 0.6940, 0.1250), \
                 linewidth = 2, label = 'Morning')
        plt.plot(cursor_time, avg_curs_sig_noon[ii], color = (0.5, 0, 0.5), \
                 linewidth = 2, label = 'Afternoon')
        
        if 'gocue' in event:
            # Solid dark green line indicating the aligned time
            plt.plot([0, 0], [y_limits[0], y_limits[1]], \
                     color = (0, 0.5, 0), linewidth = font_specs.plot_line_size)
            # Dotted dark green line indicating beginning of measured window
            plt.plot([-time_before_gocue, -time_before_gocue], \
                     [y_limits[0], y_limits[1]], color = (0, 0.5, 0), \
                     linewidth = font_specs.plot_line_size, linestyle = 'dashed')
        elif 'end' in event:
            # Solid red line indicating the aligned time
            plt.plot([0, 0], [y_limits[0], y_limits[1]], \
                     color = 'r', linewidth = font_specs.plot_line_size)
            # Dotted red line indicating beginning of measured window
            plt.plot([-time_before_end, -time_before_end], \
                     [y_limits[0], y_limits[1]], color = 'r', \
                     linewidth = font_specs.plot_line_size, linestyle = 'dashed')
        
        # Setting the axis limits
        if 'gocue' in event:
            plt.xlim([-before_event + 2, after_event])
        elif 'end' in event:
            plt.xlim([-before_event, after_event - 2])
        else:
            plt.xlim([-before_event + 1, after_event - 1])
        plt.ylim([y_limits[0], y_limits[1]])
        
        # Define the labels
        if signal_choice == 'Pos':
            signal_label = 'Wrist Position'
        elif signal_choice == 'Vel':
            signal_label = 'Wrist Velocity'
        elif signal_choice == 'Acc':
            signal_choice = 'Wrist Acceleration'
            
        # Labeling the axis
        plt.ylabel(signal_label, fontname = font_specs.font_name, \
                   fontsize = font_specs.label_font_size)
        plt.xlabel('Time (sec.)', fontname = font_specs.font_name, \
                   fontsize = font_specs.label_font_size)
            
        # Titling the plot
        title_string = 'Mean ' + signal_label + ': ' + str(target_dirs_morn[ii]) + \
            '°, TgtCenter at ' + str(target_centers_morn[ii])
        plt.title(title_string, fontname = font_specs.font_name, \
                  fontsize = font_specs.title_font_size, fontweight = 'bold')
        
        # Define the legend location
        if 'gocue' in event:
            legend_location = 'upper left'
        else:
            legend_location = 'upper right'
    
        # Legend
        legend_font = fm.FontProperties(family = font_specs.font_name, size = font_specs.legend_font_size)
        plt.legend(prop = legend_font)
        plt.legend(frameon = False, loc = legend_location)
        
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
        
    
                    
                    
                    
                    