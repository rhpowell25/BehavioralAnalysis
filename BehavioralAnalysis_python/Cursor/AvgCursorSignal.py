#%% Import basic packages
import numpy as np
import math
import matplotlib.pyplot as plt
from Plot_Specs import Font_Specs
    
#%% File Description:

# This function finds the average cursor position, velocity, or 
# acceleration per target direction / distance of all the succesful trials 
# in an xds file. It can also plot that signal.
# If you set Plot_Figs to 0, the figure will not be plotted.
# If you set Save_Figs to 0, the figure will not be saved to your desktop.
#
# -- Inputs --
# xds: the xds file
# signal_choice: 'Pos', 'Vel', or 'Acc'
# event: 'trial_gocue', 'window_trial_gocue', 'trial_end', 
# 'window_trial_end', 'force_onset', 'window_force_onset', 'force_max', 
# 'window_force_max', 'window_force_deriv', 'force_deriv', 'cursor_onset', 
# 'window_cursor_onset', 'cursor_veloc', 'window_cursor_veloc', 
# 'cursor_acc', 'window_cursor_acc', 'EMG_max', 'window_EMG_max', 
# 'task_onset', or 'window_task_onset'
# unit_name: 'elec1_1', #, or NaN
# cursor_YLims: [ymax, ymin]
# Plot_Figs: 1 or 0
# Save_Figs: 'pdf', 'png', 'fig', or 0

def AvgCursorSignal(xds, signal_choice, event, unit_name, cursor_YLims, Plot_Figs, Save_Figs):
        
    #%% Extract the target directions & centers
    from Identify_Targets import Identify_Targets
    target_vars = Identify_Targets(xds)
    target_dirs = target_vars.target_dirs
    target_centers = target_vars.target_centers
    
    #%% Basic settings, some variable extractions, & definitions
    
    # Event lengths
    before_event = 3
    after_event = 3
    
    # Window to calculate the max firing rate
    window_size = 0.1
    
    if 'window' not in event:
        max_fr_time = 0
        
    if 'gocue' in event or 'force_onset' in event:
        # Define the window for the baseline phase
        time_before_gocue = 0.4
    elif 'end' in event:
        # Define the window for the movement phase
        time_before_end = xds._lab_data__meta['TgtHold']
    
    # Extract the cursor signal of chhoice
    if signal_choice == 'Pos':
        curs_sig = xds.curs_p
    elif signal_choice == 'Vel':
        curs_sig = xds.curs_v
    elif signal_choice == 'Acc':
        curs_sig = xds.curs_a
            
    # Font & plotting specifications
    font_specs = Font_Specs()
                    
    #%% Indexes for rewarded trials in all directions
    # Counts the number of directions used
    num_dirs = len(target_dirs)
    
    # Define the output variable
    avg_curs_sig = [[] for ii in range(num_dirs)]
    
    #%% Begin the loop through all directions
    for jj in range(num_dirs):
        
        #%% Times for rewarded trials
        from EventAlignmentTimes import EventAlignmentTimes
        if event == 'trial_gocue':
            Alignment_Times = EventAlignmentTimes(xds, 'NaN', 'NaN', event)
        else:
            Alignment_Times = EventAlignmentTimes(xds, target_dirs[jj], target_centers[jj], event)
            
        if 'window' in event:
            # Run the firing rate window function
            from EventWindow import Event_Window
            Event_Window_Vars = Event_Window(xds, unit_name, target_dirs[jj], target_centers[jj], event)
            max_fr_time = Event_Window_Vars.max_fr_time
            
        #%% Extracting the cursor signal during succesful trials
        
        # Cursor signal measured during each succesful trial
        rewarded_curs_sig = [[] for ii in range(len(Alignment_Times))]
        for ii in range(len(Alignment_Times)):
            temp_start = 0
            alignment_idx = np.where(xds.time_frame == Alignment_Times[ii])[0]
            alignment_start_idx = int(alignment_idx - (before_event / xds.bin_width))
            if alignment_start_idx < 0:
                temp_start = alignment_start_idx
                alignment_start_idx = 1
            alignment_end_idx = int(alignment_idx + (after_event / xds.bin_width))
            rewarded_curs_sig[ii] = curs_sig[alignment_start_idx : alignment_end_idx]
            if temp_start != 0:
                nan_stack = np.nan * np.ones(shape = (abs(temp_start), 2))
                rewarded_curs_sig[ii] = np.vstack((nan_stack, \
                                                curs_sig[alignment_start_idx : alignment_end_idx]))
            
        #%% Find the vector sum of the cursor signal
        z_curs_sig = [[] for ii in range(len(Alignment_Times))]
        for ii in range(len(Alignment_Times)):
            z_curs_sig[ii] = np.zeros(len(rewarded_curs_sig[ii]))
            for dd in range(len(rewarded_curs_sig[ii])):
                z_curs_sig[ii][dd] = math.sqrt(rewarded_curs_sig[ii][dd][0]**2 + \
                                                 rewarded_curs_sig[ii][dd][1]**2)
          
        #%% Put all the recomposed cursor signals in a single matrix
        all_z_curs_sig = np.zeros((len(z_curs_sig[0]), len(z_curs_sig)))
        for ii in range(len(z_curs_sig)):
            all_z_curs_sig[:,ii] = z_curs_sig[ii]
            
        #%% Average the recomposed cursor signals
        avg_curs_sig[jj] = np.zeros(len(all_z_curs_sig))
        for ii in range(len(avg_curs_sig[jj])):
            avg_curs_sig[jj][ii] = np.mean(all_z_curs_sig[ii,:])
            
        #%% Find the standard dev of the recomposed cursor signals
        std_z_curs_sig = np.zeros(len(all_z_curs_sig))
        for ii in range(len(avg_curs_sig[jj])):
            std_z_curs_sig[ii] = np.std(all_z_curs_sig[ii,:])
            
        #%% Defint the absolute timing
        absolute_timing = np.linspace(-before_event, after_event, len(rewarded_curs_sig[0]))
            
        #%% Plot the average & standard deviation cursor signal
            
        if Plot_Figs == 1:

            plt.figure(figsize = (11, 4))
            
            # Average
            plt.plot(absolute_timing, avg_curs_sig[jj], 'k', linewidth = 2)
            # Standard Dev
            plt.plot(absolute_timing, avg_curs_sig[jj] + std_z_curs_sig, \
                 'r', linewidth = 1, linestyle = 'dashed')
            plt.plot(absolute_timing, avg_curs_sig[jj] - std_z_curs_sig, \
                 'r', linewidth = 1, linestyle = 'dashed')
                
            # Setting the axis limits
            if 'gocue' in event:
                plt.xlim([-before_event + 2, after_event])
            elif 'end' in event:
                plt.xlim([-before_event, after_event - 2])
            else:
                plt.xlim([-before_event + 1, after_event - 1])
            plt.ylim([cursor_YLims[1], cursor_YLims[0]])
            
            if 'gocue' in event:
                # Solid dark green line indicating the aligned time
                plt.plot([0, 0], [cursor_YLims[0], cursor_YLims[1]], \
                         color = (0, 0.5, 0), linewidth = font_specs.plot_line_size)
                # Dotted dark green line indicating beginning of measured window
                plt.plot([-time_before_gocue, -time_before_gocue], \
                         [cursor_YLims[0], cursor_YLims[1]], color = (0, 0.5, 0), \
                         linewidth = font_specs.plot_line_size, linestyle = 'dashed')
            elif 'end' in event:
                # Solid red line indicating the aligned time
                plt.plot([0, 0], [cursor_YLims[0], cursor_YLims[1]], \
                         color = 'r', linewidth = font_specs.plot_line_size)
                # Dotted red line indicating beginning of measured window
                plt.plot([-time_before_end, -time_before_end], \
                         [cursor_YLims[0], cursor_YLims[1]], color = 'r', \
                         linewidth = font_specs.plot_line_size, linestyle = 'dashed')
                    
            if 'window' in event:
                # Dotted purple line indicating beginning of measured window
                plt.plot([max_fr_time - window_size, max_fr_time - window_size], \
                         [cursor_YLims[0], cursor_YLims[1]], color = (.5, 0, .5), \
                         linewidth = font_specs.plot_line_size, linestyle = 'dashed')
                # Dotted purple line indicating end of measured window
                plt.plot([max_fr_time + window_size, max_fr_time + window_size], \
                         [cursor_YLims[0], cursor_YLims[1]], color = (.5, 0, .5), \
                         linewidth = font_specs.plot_line_size, linestyle = 'dashed')
            elif 'trial_gocue' not in event and 'trial_end' not in event:
                # Dotted red line indicating beginning of measured window
                plt.plot([-0.1, -0.1], [cursor_YLims[0], cursor_YLims[1]], color = 'r', \
                         linewidth = font_specs.plot_line_size, linestyle = 'dashed')
                plt.plot([0.1, 0.1], [cursor_YLims[0], cursor_YLims[1]], color = 'r', \
                         linewidth = font_specs.plot_line_size, linestyle = 'dashed')
            
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
            title_string = 'Mean ' + signal_label + ': ' + str(target_dirs[jj]) + \
                '°, TgtCenter at ' + str(target_centers[jj])
            plt.title(title_string, fontname = font_specs.font_name, \
                      fontsize = font_specs.title_font_size, fontweight = 'bold')
            
            # End the event after one loop if showing baseline firing rates
            if event == 'trial_gocue':
                return avg_curs_sig
            
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
            
    return avg_curs_sig
            
            
            
                    
                    
                    
                    
                    
                    
                    