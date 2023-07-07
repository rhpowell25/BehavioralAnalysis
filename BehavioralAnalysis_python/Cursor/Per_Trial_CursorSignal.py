#%% Import basic packages
import numpy as np
import math
import matplotlib.pyplot as plt
from Plot_Specs import Font_Specs
    
#%% File Description:

# This function plots the trial-by-trial cursor position, velocity, or 
# acceleration per target direction / distance of all the succesful trials 
# in an xds file. 
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
# Save_Figs: 'pdf', 'png', 'fig', or 0

def Per_Trial_CursorSignal(xds, signal_choice, event, unit_name, cursor_YLims, Save_Figs):
        
    #%% End the function if there is no Y-Limit
    
    # End the function if there is no Y-Limit
    if cursor_YLims == 'NaN':
        print("There is no Y-Limit")
        return
    
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
    
    #%% Begin the loop through all directions
    for jj in range(num_dirs):
        
        #%% Times for rewarded trials
        from EventAlignmentTimes import EventAlignmentTimes
        if event == 'trial_gocue':
            rewarded_gocue_time = EventAlignmentTimes(xds, 'NaN', 'NaN', 'trial_gocue')
            rewarded_end_time = EventAlignmentTimes(xds, 'NaN', 'NaN', 'trial_end')
            Alignment_Times = EventAlignmentTimes(xds, 'NaN', 'NaN', event)
        else:
            rewarded_gocue_time = EventAlignmentTimes(xds, target_dirs[jj], target_centers[jj], 'trial_gocue')
            rewarded_end_time = EventAlignmentTimes(xds, target_dirs[jj], target_centers[jj], 'trial_end')
            Alignment_Times = EventAlignmentTimes(xds, target_dirs[jj], target_centers[jj], event)
            
        if 'window' in event:
            # Run the firing rate window function
            from EventWindow import Event_Window
            Event_Window_Vars = Event_Window(xds, unit_name, target_dirs[jj], target_centers[jj], event)
            max_fr_time = Event_Window_Vars.max_fr_time
          
        #%% Times between events
        # Find time between the go-cue and reward
        gocue_to_event = Alignment_Times - rewarded_gocue_time
        event_to_end = rewarded_end_time - Alignment_Times
    
        #%% Extracting the cursor signal during succesful trials
        
        # Cursor signal measured during each succesful trial
        rewarded_curs_sig = [[] for ii in range(len(Alignment_Times))]
        # Time points during each succesful trial 
        timings = [[] for ii in range(len(Alignment_Times))]
        for ii in range(len(Alignment_Times)):
            temp_start = 0
            alignment_idx = np.where(xds.time_frame == Alignment_Times[ii])[0]
            alignment_start_idx = int(alignment_idx - (before_event / xds.bin_width))
            if alignment_start_idx < 0:
                temp_start = alignment_start_idx
                alignment_start_idx = 1
            alignment_end_idx = int(alignment_idx + (after_event / xds.bin_width))
            rewarded_curs_sig[ii] = curs_sig[alignment_start_idx : alignment_end_idx]
            timings[ii] = xds.time_frame[alignment_start_idx : alignment_end_idx]
            if temp_start != 0:
                nan_stack = np.nan * np.ones(shape = (abs(temp_start), 2))
                rewarded_curs_sig[ii] = np.vstack((nan_stack, \
                                                curs_sig[alignment_start_idx : alignment_end_idx]))
                timings[ii] = np.vstack((nan_stack, \
                                         xds.time_frame[alignment_start_idx : alignment_end_idx]))
            
        #%% Find the vector sum of the cursor signal
        z_curs_sig = [[] for ii in range(len(Alignment_Times))]
        for ii in range(len(Alignment_Times)):
            z_curs_sig[ii] = np.zeros(len(rewarded_curs_sig[ii]))
            for dd in range(len(rewarded_curs_sig[ii])):
                z_curs_sig[ii][dd] = math.sqrt(rewarded_curs_sig[ii][dd][0]**2 + \
                                                 rewarded_curs_sig[ii][dd][1]**2)
            
        #%% Defint the absolute timing
        absolute_timing = np.linspace(-before_event, after_event, len(rewarded_curs_sig[0]))
            
        #%% Plot the decomposed individual cursor signals on the top

        curs_sig_fig, ax = plt.subplots(figsize = (8, 8))
        ax.axis('off')
        
        dec_curs_sig = curs_sig_fig.add_subplot(2,1,1) # Top plot
        com_curs_sig = curs_sig_fig.add_subplot(2,1,2) # Bottom plot
        
        for ii in range(len(rewarded_gocue_time)):
            dec_curs_sig.plot(absolute_timing, rewarded_curs_sig[ii][:,0], \
                              'b', linewidth = .2, linestyle = 'dashed')
            dec_curs_sig.plot(absolute_timing, rewarded_curs_sig[ii][:,1], \
                              'g', linewidth = .2, linestyle = 'dashed')
            
        # Setting the axis limits
        if 'gocue' in event:
            dec_curs_sig.axis(xmin = -before_event + 2, xmax = after_event)
        elif 'end' in event:
            dec_curs_sig.axis(xmin = -before_event, xmax = after_event - 2)
        else:
            dec_curs_sig.axis(xmin = -before_event + 1, xmax = after_event - 1)
        dec_curs_sig.axis(ymin = cursor_YLims[1], ymax = cursor_YLims[0])
        
        if 'gocue' in event:
            # Solid dark green line indicating the aligned time
            dec_curs_sig.plot([0, 0], [cursor_YLims[0], cursor_YLims[1]], \
                     color = (0, 0.5, 0), linewidth = font_specs.plot_line_size)
            # Dotted dark green line indicating beginning of measured window
            dec_curs_sig.plot([-time_before_gocue, -time_before_gocue], \
                     [cursor_YLims[0], cursor_YLims[1]], color = (0, 0.5, 0), \
                     linewidth = font_specs.plot_line_size, linestyle = 'dashed')
        elif 'end' in event:
            # Solid red line indicating the aligned time
            dec_curs_sig.plot([0, 0], [cursor_YLims[0], cursor_YLims[1]], \
                     color = 'r', linewidth = font_specs.plot_line_size)
            # Dotted red line indicating beginning of measured window
            dec_curs_sig.plot([-time_before_end, -time_before_end], \
                     [cursor_YLims[0], cursor_YLims[1]], color = 'r', \
                     linewidth = font_specs.plot_line_size, linestyle = 'dashed')
                
        if 'window' in event:
            # Dotted purple line indicating beginning of measured window
            dec_curs_sig.plot([max_fr_time - window_size, max_fr_time - window_size], \
                     [cursor_YLims[0], cursor_YLims[1]], color = (.5, 0, .5), \
                     linewidth = font_specs.plot_line_size, linestyle = 'dashed')
            # Dotted purple line indicating end of measured window
            dec_curs_sig.plot([max_fr_time + window_size, max_fr_time + window_size], \
                     [cursor_YLims[0], cursor_YLims[1]], color = (.5, 0, .5), \
                     linewidth = font_specs.plot_line_size, linestyle = 'dashed')
        elif 'trial_gocue' not in event and 'trial_end' not in event:
            # Dotted red line indicating beginning of measured window
            dec_curs_sig.plot([-0.1, -0.1], [cursor_YLims[0], cursor_YLims[1]], color = 'r', \
                     linewidth = font_specs.plot_line_size, linestyle = 'dashed')
            dec_curs_sig.plot([0.1, 0.1], [cursor_YLims[0], cursor_YLims[1]], color = 'r', \
                     linewidth = font_specs.plot_line_size, linestyle = 'dashed')
        
        # Define the labels
        if signal_choice == 'Pos':
            signal_label = 'Wrist Position'
        elif signal_choice == 'Vel':
            signal_label = 'Wrist Velocity'
        elif signal_choice == 'Acc':
            signal_choice = 'Wrist Acceleration'
            
        # Labeling the axis
        dec_curs_sig.set_ylabel(signal_label, fontname = font_specs.font_name, \
                   fontsize = font_specs.label_font_size)
        dec_curs_sig.set_xlabel('Time (sec.)', fontname = font_specs.font_name, \
                   fontsize = font_specs.label_font_size)
      
        # Titling the plot
        title_string = 'Decomposed ' + signal_label + ': ' + str(target_dirs[jj]) + \
            '°, TgtCenter at ' + str(target_centers[jj])
        dec_curs_sig.set_title(title_string, fontname = font_specs.font_name, \
                  fontsize = font_specs.title_font_size, fontweight = 'bold')
        
        # Plot the individual cursor signals on the bottom
        
        for ii in range(len(rewarded_gocue_time)):
            com_curs_sig.plot(absolute_timing, z_curs_sig[ii], \
                              linewidth = .2)
                
        for ii in range(len(rewarded_gocue_time)):
            cursor_gocue_idx = timings[ii] == rewarded_gocue_time[ii]
            cursor_end_idx = timings[ii] == rewarded_end_time[ii]
            # Plot the go-cues as dark green dots
            
            
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

            
            
            
                    
                    
                    
                    
                    
                    
                    