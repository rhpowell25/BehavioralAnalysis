#%% Import basic packages
import numpy as np
import math
import scipy.stats as stats
import matplotlib.pyplot as plt
from Plot_Specs import Font_Specs

class TgtHold_CursorPos_Stats():
    
    #%% File Description:

    # This function finds changes in cursor position during the target hold 
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
        print('TtgHold Cursor Position Statistics:')
        
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
            
            rewarded_end_time_morn = EventAlignmentTimes(xds_morn, target_dirs_morn[jj], target_centers_morn[jj], 'trial_end')
            rewarded_end_time_noon = EventAlignmentTimes(xds_noon, target_dirs_noon[jj], target_centers_noon[jj], 'trial_end')
                              
            #%% Define the output variables
            Morn_TgtHold_CursorPos = np.zeros(num_dirs)
            Err_Morn_TgtHold_CursorPos = np.zeros(num_dirs)
            Noon_TgtHold_CursorPos = np.zeros(num_dirs)
            Err_Noon_TgtHold_CursorPos = np.zeros(num_dirs)
            TgtHold_CursorPos_p_values = np.zeros(num_dirs)
            TgtHold_CursorPos_perc_changes = np.zeros(num_dirs)

            #%% Cursor position aligned to trial end
            
            idx_length = int(TgtHold_time / xds_morn.bin_width)
            
            # Cursor position during each succesful trial
            cursor_p_morn  = [[] for ii in range(len(rewarded_end_time_morn))] # Cursor position at the start
            for ii in range(len(rewarded_end_time_morn)):
                idx_morn = np.argwhere(xds_morn.time_frame == rewarded_end_time_morn[ii]).reshape(-1)
                cursor_p_morn[ii] = xds_morn.curs_p[range(idx_morn[0] - idx_length, idx_morn[0])]
                
            cursor_p_noon  = [[] for ii in range(len(rewarded_end_time_noon))] # Cursor position at the start
            for ii in range(len(rewarded_end_time_noon)):
                idx_noon = np.argwhere(xds_noon.time_frame == rewarded_end_time_noon[ii]).reshape(-1)
                cursor_p_noon[ii] = xds_noon.curs_p[range(idx_noon[0] - idx_length, idx_noon[0])]
                
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
                
            #%% Putting all succesful trials in one array
            all_trials_z_cursor_p_morn = np.zeros((len(z_cursor_p_morn[0]), len(rewarded_end_time_morn)))
            for ii in range(len(rewarded_end_time_morn)):
                all_trials_z_cursor_p_morn[:,ii] = z_cursor_p_morn[ii]
                  
            all_trials_z_cursor_p_noon = np.zeros((len(z_cursor_p_noon[0]), len(rewarded_end_time_noon)))
            for ii in range(len(rewarded_end_time_noon)):
                all_trials_z_cursor_p_noon[:,ii] = z_cursor_p_noon[ii]
                    
            #%% Normalizing the average cursor position
            from Multi_Session_NormalizeCursor import Multi_Session_NormalizeCursor
            Cursor_Norm_Factor = Multi_Session_NormalizeCursor(xds_morn, xds_noon, norm_cursor)
            
            all_trials_z_cursor_p_morn = all_trials_z_cursor_p_morn / Cursor_Norm_Factor * 100
            all_trials_z_cursor_p_noon = all_trials_z_cursor_p_noon / Cursor_Norm_Factor * 100
            
            #%% Calculating average cursor position (Average per trial)
            per_trial_avg_z_cursor_p_morn = np.zeros(len(cursor_p_morn))
            per_trial_avg_z_cursor_p_noon = np.zeros(len(cursor_p_noon))
            for ii in range(len(cursor_p_morn)):
                per_trial_avg_z_cursor_p_morn[ii] = np.mean(all_trials_z_cursor_p_morn[:,ii])
            for ii in range(len(cursor_p_noon)):
                per_trial_avg_z_cursor_p_noon[ii] = np.mean(all_trials_z_cursor_p_noon[:,ii])
                    
            #%% Do the statistics on the cursor position
            TgtHold_CursorPos_p_values[jj] = \
                stats.ttest_ind(per_trial_avg_z_cursor_p_morn, per_trial_avg_z_cursor_p_noon)[1]
            # Average cursor position
            Morn_TgtHold_CursorPos[jj] = np.mean(per_trial_avg_z_cursor_p_morn)
            Noon_TgtHold_CursorPos[jj] = np.mean(per_trial_avg_z_cursor_p_noon)
            # Standard deviation
            TgtHold_cursor_p_morn_std = np.std(per_trial_avg_z_cursor_p_morn)
            TgtHold_cursor_p_noon_std = np.std(per_trial_avg_z_cursor_p_noon)
            # Standard error
            Err_Morn_TgtHold_CursorPos[jj] = \
                TgtHold_cursor_p_morn_std / math.sqrt(len(per_trial_avg_z_cursor_p_morn))
            Err_Noon_TgtHold_CursorPos[jj] = \
                TgtHold_cursor_p_noon_std / math.sqrt(len(per_trial_avg_z_cursor_p_noon))
            # Cursor position percent change
            if Cursor_Norm_Factor == 1:
                TgtHold_CursorPos_perc_changes[jj] = \
                    (Noon_TgtHold_CursorPos[jj] - Morn_TgtHold_CursorPos[jj]) / 100
            else:
                TgtHold_CursorPos_perc_changes[jj] = \
                    (Noon_TgtHold_CursorPos[jj] - Morn_TgtHold_CursorPos[jj]) / abs(Morn_TgtHold_CursorPos[jj])
            
            #%% Save the necessary variables
            self.TgtHold_CursorPos_p_values = TgtHold_CursorPos_p_values
            self.TgtHold_CursorPos_perc_changes = TgtHold_CursorPos_perc_changes
            self.Morn_TgtHold_CursorPos = Morn_TgtHold_CursorPos
            self.Err_Morn_TgtHold_CursorPos = Err_Morn_TgtHold_CursorPos
            self.Noon_TgtHold_CursorPos = Noon_TgtHold_CursorPos
            self.Err_Noon_TgtHold_CursorPos = Err_Noon_TgtHold_CursorPos
 
        #%% Plot the cursor position box & whisker plots
        if Plot_Figs == 1:
    
            # Combining the morning & afternoon into one array
            avg_cursor_p = [[] for ii in range(2)]
            avg_cursor_p[0] = per_trial_avg_z_cursor_p_morn
            avg_cursor_p[1] = per_trial_avg_z_cursor_p_noon
            
            TgtHold_cursor_p_mean = np.hstack((Morn_TgtHold_CursorPos, Noon_TgtHold_CursorPos))
            TgtHold_cursor_p_std = np.hstack((TgtHold_cursor_p_morn_std, TgtHold_cursor_p_noon_std))
        
            # Finding the min and max for the y-axis limit
            TgtHold_y_min = TgtHold_cursor_p_mean - TgtHold_cursor_p_std
            TgtHold_y_max = TgtHold_cursor_p_mean + TgtHold_cursor_p_std
            
            # Define the figure
            cursor_pos_fig = plt.figure(figsize=(8, 6))
            box_plot = cursor_pos_fig.add_subplot(2,1,1) # Top plot
            morn_cursor_pos = cursor_pos_fig.add_subplot(2,2,3) # Bottom left plot
            noon_cursor_pos = cursor_pos_fig.add_subplot(2,2,4)
            
            # Boxplot
            box_plot.boxplot(avg_cursor_p)
            box_plot.set_xticklabels('')
            
            # Setting the y-axis limits
            y_max = max(TgtHold_y_max)
            y_min = min(TgtHold_y_min)
            box_plot.axis(ymin = y_min - 1.5*abs(y_min), ymax = y_max + 1.5*abs(y_max))
            # Setting the x-axis limits
            box_plot.axis(xmin = 0.5, xmax = 2.5)
            
            # Titling the plot
            title_string = 'TgtHold Wrist Position, ' + str(round(target_dirs_morn[jj])) \
                + '°, TgtCenter at ' + str(round(target_centers_morn[jj]))
            box_plot.set_title(title_string, fontname = font_specs.font_name, 
                         fontsize = font_specs.title_font_size + 5, fontweight = 'bold')
                    
            # Get the top subplot title for saving
            if Save_Figs != 0:
                save_title[jj] = title_string
            
            # Annotation of the p-value
            if round(TgtHold_CursorPos_p_values[jj], 3) > 0:
                box_plot.text(0.1, 0.85, 'p = ' + str(round(TgtHold_CursorPos_p_values[jj], 3)), \
                         verticalalignment = 'center', horizontalalignment = 'center', \
                             transform = box_plot.transAxes, fontname = font_specs.font_name, \
                                 fontsize = font_specs.legend_font_size)
            if round(TgtHold_CursorPos_p_values[jj], 3) == 0:
                box_plot.text(0.1, 0.85, 'p < 0.001', verticalalignment = 'center', \
                         horizontalalignment = 'center', transform = box_plot.transAxes, \
                             fontname = font_specs.font_name, fontsize = font_specs.legend_font_size)
            
            # Annotation of the percent change
            if round(TgtHold_CursorPos_perc_changes[jj], 3) != 0:
                box_plot.text(0.875, 0.85, 'Δ% = ' + str(round(TgtHold_CursorPos_perc_changes[jj], 3)), \
                         verticalalignment = 'center', horizontalalignment = 'center', \
                             transform = box_plot.transAxes, fontname = font_specs.font_name, \
                                 fontsize = font_specs.legend_font_size)
            if round(TgtHold_CursorPos_perc_changes[jj], 3) == 0:
                box_plot.text(0.875, 0.85, 'Δ% ≈ 0', \
                         verticalalignment = 'center', horizontalalignment = 'center', \
                             transform = box_plot.transAxes, fontname = font_specs.font_name, \
                                 fontsize = font_specs.legend_font_size)
            
        # Bottom plots
        per_trial_ymax_morn = np.amax(all_trials_z_cursor_p_morn)
        per_trial_ymax_noon = np.amax(all_trials_z_cursor_p_noon)
        per_trial_ymin_morn = np.amin(all_trials_z_cursor_p_morn)
        per_trial_ymin_noon = np.amin(all_trials_z_cursor_p_noon)
        y_max = max(per_trial_ymax_morn, per_trial_ymax_noon)
        y_min = min(per_trial_ymin_morn, per_trial_ymin_noon)
        
        # Set the title
        morn_cursor_pos.set_title('Morning', fontname = font_specs.font_name, 
                       fontsize = font_specs.title_font_size + 5, fontweight = 'bold')
        morn_cursor_pos.plot(all_trials_z_cursor_p_morn)
        
        # Annotation of the number of trials
        morn_cursor_pos.text(0.2, 0.85, 'n = ' + str(len(rewarded_end_time_morn)), \
                             verticalalignment = 'center', horizontalalignment = 'center', \
                                 transform = morn_cursor_pos.transAxes, fontname = font_specs.font_name, \
                                     fontsize = font_specs.legend_font_size)
        
        # Setting the axis limits
        morn_cursor_pos.axis(ymin = y_min - abs(y_min/8), ymax = y_max + abs(y_max/8))
        morn_cursor_pos.axis(xmin = 0, xmax = idx_length)
        
        # Set the title
        noon_cursor_pos.set_title('Afternoon', fontname = font_specs.font_name, 
                       fontsize = font_specs.title_font_size + 5, fontweight = 'bold')
        noon_cursor_pos.plot(all_trials_z_cursor_p_noon)
        
        # Annotation of the number of trials
        noon_cursor_pos.text(1.35, 0.85, 'n = ' + str(len(rewarded_end_time_noon)), \
                             verticalalignment = 'center', horizontalalignment = 'center', \
                                 transform = morn_cursor_pos.transAxes, fontname = font_specs.font_name, \
                                     fontsize = font_specs.legend_font_size)
        
        # Setting the axis limits
        noon_cursor_pos.axis(ymin = y_min - abs(y_min/8), ymax = y_max + abs(y_max/8))
        noon_cursor_pos.axis(xmin = 0, xmax = idx_length)
        
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
        
        
        
        
        
        
        
        
        
        
        
        
        
            
                    