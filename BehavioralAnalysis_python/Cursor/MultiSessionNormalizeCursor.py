# -*- coding: utf-8 -*-
import numpy as np
import math

class MultiSession_NormalizeCursor():
    
    def Multi_Session_Normalize_Cursor(xds_morn, xds_noon, muscle_groups, norm_perctile, norm_cursor):
        
        #%% End the function if you are not normalizing the EMG
        if norm_cursor == 0:
            print('Cursor Position Will Not Be Normalized')
            Cursor_Norm_Factor = 1
            return Cursor_Norm_Factor
        else:
            print('Cursor Will Be Normalized to the ' + str(norm_perctile) + 'th percentile')

        #%% Concatenate all the morning & afternoon information
        
        # Time frame
        noon_time_frame = xds_noon.time_frame + xds_morn.time_frame[-1]
        time_frame = np.concatenate((xds_morn.time_frame, noon_time_frame), axis = 0)
        
        # Bin width
        bin_width = xds_morn.bin_width
        
        # Cursor
        cat_Cursor = np.concatenate((xds_morn.curs_p, xds_noon.curs_p), axis = 0)
        
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
            
        #%% Pull the cursor position corresponding to the extracted time frames
        
        # Pull all the cursor position
        cursor_p = [[] for ii in range(len(total_rewarded_idx))]
        for ii in range(len(total_rewarded_idx)):
            if rewarded_end_idx[ii] + (2/bin_width) > len(cat_Cursor):
                cursor_p[ii] = cat_Cursor[int(rewarded_start_idx[ii]) : -1, :]
            else:
                cursor_p[ii] = cat_Cursor[int(rewarded_start_idx[ii]) : int(rewarded_start_idx[ii] + (2/bin_width)), :]
         
        #%% Recompose the cursor position
        
        z_Cursor = [[] for ii in range(len(total_rewarded_idx))]
        # Loops through cursor position
        for ii in range(len(total_rewarded_idx)):
            z_Cursor[ii] = np.zeros(len(cursor_p[ii]))
            for dd in range(len(z_Cursor[ii])):
                z_Cursor[ii][dd] = math.sqrt(cursor_p[ii][dd][0]**2 + cursor_p[ii][dd][1]**2)
                        
        #%% Finding the max cursor position per trial
        
        max_Cursor_pertrial = np.zeros(len(total_rewarded_idx))
        for ii in range(len(total_rewarded_idx)):
            max_Cursor_pertrial[ii] = max(z_Cursor[ii])
                
        #%% Find the 95th percentile of the max's
        Cursor_Norm_Factor = np.percentile(max_Cursor_pertrial, norm_perctile)
    
    
        return Cursor_Norm_Factor
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    