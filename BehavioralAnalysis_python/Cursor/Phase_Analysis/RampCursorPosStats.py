#%% Import basic packages
import numpy as np
import math
import scipy.stats as stats
import matplotlib.pyplot as plt
from Plot_Specs import Font_Specs

class Ramp_CursorPos_Stats():
    
    #%% File Description:

    # This function finds changes in cursor position during the ramp phase
    # between two XDS files & plots a boxplot & individual trial traces.
    # If you set Plot_Figs to 0, the figure will not be plotted.
    # If you set Save_Figs to 0, the figure will not be saved to your desktop.
    #
    # -- Inputs --
    # xds_morn: the first xds file
    # xds_noon: the second xds file
    # Cursor_Norm_Factor: the factor to normalize the cursor position
    # Plot_Figs: 0 or 1
    # Save_Figs: 'pdf', 'png', or 0

    def __init__(self, xds_morn, xds_noon, norm_cursor, Plot_Figs, Save_Figs):
        
        #%% Display the function being used
        print('Ramp Cursor Position Statistics:')
        
        #%% Basic Settings, some variable extractions, & definitions
        
        # Font & plotting specifications
        font_specs = Font_Specs()
        
        # Define the window for the baseline phase
        TgtHold_time = xds_morn._lab_data__meta['TgtHold']

        #%% Extract the target directions & centers
        from Identify_Targets import Identify_Targets
        target_vars_morn = Identify_Targets(xds_morn)
        target_dirs_morn = target_vars_morn.target_dirs
        target_centers_morn = target_vars_morn.target_centers
        target_vars_noon = Identify_Targets(xds_noon)
        target_dirs_noon = target_vars_noon.target_dirs
        target_centers_noon = target_vars_noon.target_centers
        
        #%% Check to see if both sessions use a consistent number of targets
        
        # Find matching targets between the two seesions
        from Match_Targets import Match_Targets
        Matching_Idxs = Match_Targets(target_dirs_morn, target_dirs_noon, target_centers_morn, target_centers_noon)
        Matching_Idxs_Morn = Matching_Idxs.Matching_Idxs_Morn
        Matching_Idxs_Noon = Matching_Idxs.Matching_Idxs_Noon
        
        # Only use the info of target centers conserved between morn & noon
        if not all(ii is True for ii in Matching_Idxs_Morn) or not all(ii is True for ii in Matching_Idxs_Noon):
            print('Uneven Targets Between Morning & Afternoon')
            target_centers_morn = target_centers_morn[Matching_Idxs_Morn]
            target_centers_noon = target_centers_noon[Matching_Idxs_Noon]
            target_dirs_morn = target_dirs_morn[Matching_Idxs_Morn]
            target_dirs_noon = target_dirs_noon[Matching_Idxs_Noon]
        
        #%% Settings to loop through every target direction
        
        # Counts the number of directions used
        num_dirs = len(target_dirs_morn)
        
        if Save_Figs != 0:
            save_title = [[] for ii in range(num_dirs)]
            
        #%% Begin the loop through all directions
        for jj in range(num_dirs):
            
            #%% Times for rewarded trials
            from EventAlignmentTimes import EventAlignmentTimes
            
            rewarded_onset_time_morn = EventAlignmentTimes(xds_morn, target_dirs_morn[jj], target_centers_morn[jj], 'task_onset')
            rewarded_end_time_morn = EventAlignmentTimes(xds_morn, target_dirs_morn[jj], target_centers_morn[jj], 'trial_end')
            rewarded_end_time_morn = rewarded_end_time_morn - TgtHold_time
            # Round to match the neural bin size
            rewarded_end_time_morn = np.round(rewarded_end_time_morn, abs(math.floor(math.log10(xds_morn.bin_width))))
            
            rewarded_onset_time_noon = EventAlignmentTimes(xds_noon, target_dirs_noon[jj], target_centers_noon[jj], 'task_onset')
            rewarded_end_time_noon = EventAlignmentTimes(xds_noon, target_dirs_noon[jj], target_centers_noon[jj], 'trial_end')
            rewarded_end_time_noon = rewarded_end_time_noon - TgtHold_time
            # Round to match the neural bin size
            rewarded_end_time_noon = np.round(rewarded_end_time_noon, abs(math.floor(math.log10(xds_noon.bin_width))))
                              
            #%% Define the output variables
            Morn_Ramp_CursorPos = np.zeros(num_dirs)
            STD_Morn_TgtHold_cursor_p = np.zeros(num_dirs)
            Err_Morn_Ramp_CursorPos = np.zeros(num_dirs)
            Noon_Ramp_CursorPos = np.zeros(num_dirs)
            STD_Noon_TgtHold_cursor_p = np.zeros(num_dirs)
            Err_Noon_Ramp_CursorPos = np.zeros(num_dirs)
            Ramp_CursorPos_p_values = np.zeros(num_dirs)
            Ramp_CursorPos_perc_changes = np.zeros(num_dirs)
            
            #%% Cursor position aligned to specified event
            # Find the rewarded times in the whole trial time frame
            rewarded_start_idx_morn = np.zeros(len(rewarded_onset_time_morn))
            for ii in range(len(rewarded_onset_time_morn)):
                rewarded_start_idx_morn[ii] = np.argwhere(xds_morn.time_frame == rewarded_onset_time_morn[ii])
                
            rewarded_end_idx_morn = np.zeros(len(rewarded_end_time_morn))
            for ii in range(len(rewarded_end_time_morn)):
                rewarded_end_idx_morn[ii] = np.argwhere(xds_morn.time_frame == rewarded_end_time_morn[ii])
            
            rewarded_start_idx_noon = np.zeros(len(rewarded_onset_time_noon))
            for ii in range(len(rewarded_onset_time_noon)):
                rewarded_start_idx_noon[ii] = np.argwhere(xds_noon.time_frame == rewarded_onset_time_noon[ii])
                
            rewarded_end_idx_noon = np.zeros(len(rewarded_end_time_noon))
            for ii in range(len(rewarded_end_time_noon)):
                rewarded_end_idx_noon[ii] = np.argwhere(xds_noon.time_frame == rewarded_end_time_noon[ii])
                
            # Cursor position during each succesful trial
            cursor_p_morn  = [[] for ii in range(len(rewarded_end_time_morn))] # Cursor position at the start
            for ii in range(len(rewarded_end_time_morn)):
                cursor_p_morn[ii] = xds_morn.curs_p[range(int(rewarded_start_idx_morn[ii]), int(rewarded_end_idx_morn[ii]))]
                
            cursor_p_noon  = [[] for ii in range(len(rewarded_end_time_noon))] # Cursor position at the start
            for ii in range(len(rewarded_end_time_noon)):
                cursor_p_noon[ii] = xds_noon.curs_p[range(int(rewarded_start_idx_noon[ii]), int(rewarded_end_idx_noon[ii]))]
                
            #%% Find the vector sum of the cursor position
            
            z_cursor_p_morn = [[] for ii in range(len(rewarded_end_time_morn))]
            for ii in range(len(rewarded_end_time_morn)):
                z_cursor_p_morn[ii] = np.zeros(len(cursor_p_morn[ii]))
                for dd in range(len(z_cursor_p_morn[ii])):
                    z_cursor_p_morn[ii][dd] = math.sqrt(cursor_p_morn[ii][dd][0]**2 + cursor_p_morn[ii][dd][1]**2)
                    
            z_cursor_p_noon = [[] for ii in range(len(rewarded_end_time_noon))]
            for ii in range(len(rewarded_end_time_noon)):
                z_cursor_p_noon[ii] = np.zeros(len(cursor_p_noon[ii]))
                for dd in range(len(z_cursor_p_noon[ii])):
                    z_cursor_p_noon[ii][dd] = math.sqrt(cursor_p_noon[ii][dd][0]**2 + cursor_p_noon[ii][dd][1]**2)
                    
            #%% Normalizing the average cursor position
            from Multi_Session_NormalizeCursor import Multi_Session_NormalizeCursor
            Cursor_Norm_Factor = Multi_Session_NormalizeCursor(xds_morn, xds_noon, norm_cursor)
            
            for ii in range(len(rewarded_onset_time_morn)):
                z_cursor_p_morn[ii] = z_cursor_p_morn[ii] / Cursor_Norm_Factor * 100
                
            for ii in range(len(rewarded_onset_time_noon)):
                z_cursor_p_noon[ii] = z_cursor_p_noon[ii] / Cursor_Norm_Factor * 100
            
            #%% Calculating average cursor position (Average per trial)
            per_trial_avg_z_cursor_p_morn = np.zeros(len(cursor_p_morn))
            per_trial_avg_z_cursor_p_noon = np.zeros(len(cursor_p_noon))
            for ii in range(len(cursor_p_morn)):
                per_trial_avg_z_cursor_p_morn[ii] = np.nanmean(z_cursor_p_morn[ii])
            for ii in range(len(cursor_p_noon)):
                per_trial_avg_z_cursor_p_noon[ii] = np.nanmean(z_cursor_p_noon[ii])
                    
            #%% Do the statistics on the cursor position
            Ramp_CursorPos_p_values[jj] = \
                stats.ttest_ind(per_trial_avg_z_cursor_p_morn, per_trial_avg_z_cursor_p_noon, nan_policy = 'omit')[1]
            # Average cursor position
            Morn_Ramp_CursorPos[jj] = np.nanmean(per_trial_avg_z_cursor_p_morn)
            Noon_Ramp_CursorPos[jj] = np.nanmean(per_trial_avg_z_cursor_p_noon)
            # Standard deviation
            STD_Morn_TgtHold_cursor_p[jj] = np.nanstd(per_trial_avg_z_cursor_p_morn)
            STD_Noon_TgtHold_cursor_p[jj] = np.nanstd(per_trial_avg_z_cursor_p_noon)
            # Standard error
            Err_Morn_Ramp_CursorPos[jj] = \
                STD_Morn_TgtHold_cursor_p[jj] / math.sqrt(len(per_trial_avg_z_cursor_p_morn))
            Err_Noon_Ramp_CursorPos[jj] = \
                STD_Noon_TgtHold_cursor_p[jj] / math.sqrt(len(per_trial_avg_z_cursor_p_noon))
            # Cursor position percent change
            if Cursor_Norm_Factor == 1:
                Ramp_CursorPos_perc_changes[jj] = \
                    (Noon_Ramp_CursorPos[jj] - Morn_Ramp_CursorPos[jj]) / 100
            else:
                Ramp_CursorPos_perc_changes[jj] = \
                    (Noon_Ramp_CursorPos[jj] - Morn_Ramp_CursorPos[jj]) / abs(Morn_Ramp_CursorPos[jj])
            
            #%% Save the necessary variables
            self.Ramp_CursorPos_p_values = Ramp_CursorPos_p_values
            self.Ramp_CursorPos_perc_changes = Ramp_CursorPos_perc_changes
            self.Morn_Ramp_CursorPos = Morn_Ramp_CursorPos
            self.Err_Morn_Ramp_CursorPos = Err_Morn_Ramp_CursorPos
            self.Noon_Ramp_CursorPos = Noon_Ramp_CursorPos
            self.Err_Noon_Ramp_CursorPos = Err_Noon_Ramp_CursorPos
 
        #%% Plot the cursor position box & whisker plots
        if Plot_Figs == 1:
    
            # Combining the morning & afternoon into one array
            avg_cursor_p = [[] for ii in range(2)]
            avg_cursor_p[0] = per_trial_avg_z_cursor_p_morn[~np.isnan(per_trial_avg_z_cursor_p_morn)]
            avg_cursor_p[1] = per_trial_avg_z_cursor_p_noon[~np.isnan(per_trial_avg_z_cursor_p_noon)]
            
            Ramp_cursor_p_mean = np.hstack((Morn_Ramp_CursorPos, Noon_Ramp_CursorPos))
            Ramp_cursor_p_std = np.hstack((STD_Morn_TgtHold_cursor_p[jj], STD_Noon_TgtHold_cursor_p[jj]))
        
            # Finding the min and max for the y-axis limit
            Ramp_y_min = Ramp_cursor_p_mean - Ramp_cursor_p_std
            Ramp_y_max = Ramp_cursor_p_mean + Ramp_cursor_p_std
            
            # Define the figure
            cursor_pos_fig = plt.figure(figsize=(8, 6))
            box_plot = cursor_pos_fig.add_subplot(2,1,1) # Top plot
            morn_cursor_pos = cursor_pos_fig.add_subplot(2,2,3) # Bottom left plot
            noon_cursor_pos = cursor_pos_fig.add_subplot(2,2,4)
            
            # Boxplot
            box_plot.boxplot(avg_cursor_p)
            box_plot.set_xticklabels('')
            
            # Setting the y-axis limits
            y_max = max(Ramp_y_max)
            y_min = min(Ramp_y_min)
            box_plot.axis(ymin = y_min - 1.5*abs(y_min), ymax = y_max + 1.5*abs(y_max))
            # Setting the x-axis limits
            box_plot.axis(xmin = 0.5, xmax = 2.5)
            
            # Titling the plot
            title_string = 'Ramp Wrist Position, ' + str(round(target_dirs_morn[jj])) \
                + '°, TgtCenter at ' + str(round(target_centers_morn[jj]))
            box_plot.set_title(title_string, fontname = font_specs.font_name, 
                         fontsize = font_specs.title_font_size + 5, fontweight = 'bold')
                    
            # Get the top subplot title for saving
            if Save_Figs != 0:
                save_title[jj] = title_string
            
            # Annotation of the p-value
            if round(Ramp_CursorPos_p_values[jj], 3) > 0:
                box_plot.text(0.1, 0.85, 'p = ' + str(round(Ramp_CursorPos_p_values[jj], 3)), \
                         verticalalignment = 'center', horizontalalignment = 'center', \
                             transform = box_plot.transAxes, fontname = font_specs.font_name, \
                                 fontsize = font_specs.legend_font_size)
            if round(Ramp_CursorPos_p_values[jj], 3) == 0:
                box_plot.text(0.1, 0.85, 'p < 0.001', verticalalignment = 'center', \
                         horizontalalignment = 'center', transform = box_plot.transAxes, \
                             fontname = font_specs.font_name, fontsize = font_specs.legend_font_size)
            
            # Annotation of the percent change
            if round(Ramp_CursorPos_perc_changes[jj], 3) != 0:
                box_plot.text(0.875, 0.85, 'Δ% = ' + str(round(Ramp_CursorPos_perc_changes[jj], 3)), \
                         verticalalignment = 'center', horizontalalignment = 'center', \
                             transform = box_plot.transAxes, fontname = font_specs.font_name, \
                                 fontsize = font_specs.legend_font_size)
            if round(Ramp_CursorPos_perc_changes[jj], 3) == 0:
                box_plot.text(0.875, 0.85, 'Δ% ≈ 0', \
                         verticalalignment = 'center', horizontalalignment = 'center', \
                             transform = box_plot.transAxes, fontname = font_specs.font_name, \
                                 fontsize = font_specs.legend_font_size)
        
        # Set the title
        morn_cursor_pos.set_title('Morning', fontname = font_specs.font_name, 
                       fontsize = font_specs.title_font_size + 5, fontweight = 'bold')
        # Plotting
        for pp in range(len(z_cursor_p_morn)):
            morn_cursor_pos.plot(z_cursor_p_morn[pp])
        
        # Annotation of the number of trials
        morn_cursor_pos.text(0.2, 0.85, 'n = ' + str(len(rewarded_end_time_morn)), \
                             verticalalignment = 'center', horizontalalignment = 'center', \
                                 transform = morn_cursor_pos.transAxes, fontname = font_specs.font_name, \
                                     fontsize = font_specs.legend_font_size)
        
        # Set the title
        noon_cursor_pos.set_title('Afternoon', fontname = font_specs.font_name, 
                       fontsize = font_specs.title_font_size + 5, fontweight = 'bold')
        # Plotting
        for pp in range(len(z_cursor_p_noon)):
            noon_cursor_pos.plot(z_cursor_p_noon[pp])
        
        # Annotation of the number of trials
        noon_cursor_pos.text(1.35, 0.85, 'n = ' + str(len(rewarded_end_time_noon)), \
                             verticalalignment = 'center', horizontalalignment = 'center', \
                                 transform = morn_cursor_pos.transAxes, fontname = font_specs.font_name, \
                                     fontsize = font_specs.legend_font_size)
        
        # Figure Saving
        if Save_Figs != 0:
            save_dir = 'C:/Users/rhpow/Desktop/'
            fig_title = save_title[jj]
            fig_title = str.replace(fig_title, ':', '')
            fig_title = str.replace(fig_title, 'vs.', 'vs')
            fig_title = str.replace(fig_title, 'mg.', 'mg')
            fig_title = str.replace(fig_title, 'kg.', 'kg')
            fig_title = str.replace(fig_title, '.', '_')
            fig_title = str.replace(fig_title, '/', '_')
            plt.savefig(save_dir + fig_title + '.' + Save_Figs)
            plt.close()
        
        
        
        
        
        
        
        
        
        
        
        
        
            
                    