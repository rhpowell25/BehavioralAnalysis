#%% Import basic packages
import numpy as np
import math
    
#%% File Description:

# This function finds the nth percentile of cursor position, velocity, or 
# acceleration between two XDS files for the purpose of normalization.
# Depending on the method, the nth percentile will be calculated using
# the entirety of both files or using only the succesful trials.
# If you set norm_cursor to 0, the Cursor_Norm_Factor will be 1
#
# -- Inputs --
# xds_morn: the first xds file
# xds_noon: the second xds file
# signal_choice: 'Pos', 'Vel', or 'Acc'
# norm_cursor: 1 or 0

def Multi_Session_NormalizeCursor(xds_morn, xds_noon, signal_choice, norm_cursor):
    
    #%% End the function if you are not normalizing the EMG
    if norm_cursor == 0:
        print('Cursor Signal Will Not Be Normalized')
        Cursor_Norm_Factor = 1
        return Cursor_Norm_Factor
    else:
        norm_prctile = 95
        print('Cursor Signal Will Be Normalized to the ' + str(norm_prctile) + 'th percentile')

    #%% What part of the cursor signal do you want to take the percentile of
    # All the cursor data ('All_Cursor')
    # The cursor in each succesful trial ('Trial_Cursor')
    norm_method = 'Trial_Cursor'
    
    #%% Basic Settings, some variable extractions, & definitions

    # Extract the cursor signal of chhoice
    if signal_choice == 'Pos':
        curs_morn = xds_morn.curs_p
        curs_noon = xds_noon.curs_p
    elif signal_choice == 'Vel':
        curs_morn = xds_morn.curs_v
        curs_noon = xds_noon.curs_v
    elif signal_choice == 'Acc':
        curs_morn = xds_morn.curs_a
        curs_noon = xds_noon.curs_a
    
    #%% Concatenate the cursor signal
    
    cat_curs = np.concatenate((curs_morn, curs_noon), axis = 0)
    
    #%% Find the vector sum of the cursor signal
    z_cat_curs = np.zeros(len(cat_curs))
    for ii in range(len(z_cat_curs)):
        z_cat_curs[ii] = math.sqrt(cat_curs[ii,0]**2 + cat_curs[ii,1]**2)
        
    #%% Finding the max cursor signal throughout both files
    
    if norm_method == 'All_Cursor':
        
        # Find the nth percentile of the cursor signal
        Cursor_Norm_Factor = np.percentile(z_cat_curs, norm_prctile)

    #%% Concatenate all the morning & afternoon information
    if norm_method == 'Trial_Cursor':
        
        # Time frame
        noon_time_frame = xds_noon.time_frame + xds_morn.time_frame[-1]
        time_frame = np.concatenate((xds_morn.time_frame, noon_time_frame), axis = 0)
        
        # Bin width
        bin_width = xds_morn.bin_width
        
        # Trial results
        trial_results = np.concatenate((xds_morn.trial_result, xds_noon.trial_result), axis = 0)
        # Trial start times
        noon_trial_start_time = xds_noon.trial_start_time + xds_morn.time_frame[-1]
        trial_start_time = np.concatenate((xds_morn.trial_start_time, noon_trial_start_time), axis = 0)
        # Trial go cue times
        noon_trial_gocue_time = xds_noon.trial_gocue_time + xds_morn.time_frame[-1]
        trial_gocue_time = np.concatenate((xds_morn.trial_gocue_time, noon_trial_gocue_time), axis = 0)
        # Trial end times
        noon_trial_end_time = xds_noon.trial_end_time + xds_morn.time_frame[-1]
        trial_end_time = np.concatenate((xds_morn.trial_end_time, noon_trial_end_time), axis = 0)
                          
        #%% Index for rewarded trials
        
        total_rewarded_idx = np.argwhere(trial_results == 'R').reshape(-1,)
        
        #%% Loops to extract only rewarded trials
        
        # Rewarded start times
        rewarded_start_time = np.zeros(len(total_rewarded_idx))
        for ii in range(len(total_rewarded_idx)):
            rewarded_start_time[ii] = trial_start_time[total_rewarded_idx[ii]]
            
        # Rewarded go-cue's
        rewarded_gocue_time = np.zeros(len(total_rewarded_idx))
        for ii in range(len(total_rewarded_idx)):
            rewarded_gocue_time[ii] = trial_gocue_time[total_rewarded_idx[ii]]
                
        # Rewarded end times
        rewarded_end_time = np.zeros(len(total_rewarded_idx))
        for ii in range(len(total_rewarded_idx)):
            rewarded_end_time[ii] = trial_end_time[total_rewarded_idx[ii]]
            
        #%% Pulling the timeframe of the succesful trials
        # Find the rewarded start times in the whole trial time frame
        rewarded_start_idx = np.zeros(len(total_rewarded_idx))
        for ii in range(len(total_rewarded_idx)):
            rewarded_start_idx[ii] = np.asarray(np.where(time_frame == rewarded_start_time[ii]))
            
        # Find the rewarded end times in the whole trial time frame
        rewarded_end_idx = np.zeros(len(total_rewarded_idx))
        for ii in range(len(total_rewarded_idx)):
            rewarded_end_idx[ii] = np.asarray(np.where(time_frame == rewarded_end_time[ii]))  
            
        timings = [[] for ii in range(len(total_rewarded_idx))]
        for ii in range(len(total_rewarded_idx)):
            timings[ii] = time_frame[int(rewarded_start_idx[ii]) : 
                 int(rewarded_end_idx[ii] + (2/bin_width))]
            
        #%% Pull the cursor signal corresponding to the extracted time frames
        
        # Cursor signal measured during each successful trial
        rewarded_curs_sig = [[] for ii in range(len(total_rewarded_idx))]
        for ii in range(len(total_rewarded_idx)):
            if rewarded_end_idx[ii] + (2/bin_width) > len(z_cat_curs):
                rewarded_curs_sig[ii] = z_cat_curs[int(rewarded_start_idx[ii]) : -1]
            else:
                rewarded_curs_sig[ii] = z_cat_curs[int(rewarded_start_idx[ii]) : int(rewarded_start_idx[ii] + (2/bin_width))]
                
        #%% Finding the max cursor signal per trial
        
        max_curs_pertrial = np.zeros(len(total_rewarded_idx))
        for ii in range(len(total_rewarded_idx)):
            max_curs_pertrial[ii] = max(rewarded_curs_sig[ii])
                
        #%% Find the 95th percentile of the max's
        Cursor_Norm_Factor = np.percentile(max_curs_pertrial, norm_prctile)


    return Cursor_Norm_Factor
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    