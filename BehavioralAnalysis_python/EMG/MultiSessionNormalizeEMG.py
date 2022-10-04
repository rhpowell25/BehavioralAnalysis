# -*- coding: utf-8 -*-
from EMGIndex import EMG_Index
import numpy as np

class MultiSession_NormalizeEMG():
    
    def Multi_Session_Normalize_EMG(xds_morn, xds_noon, muscle_groups, norm_perctile, norm_EMG):
        
        #%% Find the EMG index
        
        M = EMG_Index.EMG_Indexing(xds_morn, muscle_groups)
        
        #%% End the function if you are not normalizing the EMG
        if norm_EMG == 0:
            print('EMG Will Not Be Normalized')
            EMG_Norm_Factor = np.ones([1,len(M)])
            return EMG_Norm_Factor
        else:
            print('EMG Will Be Normalized to the ' + str(norm_perctile) + 'th percentile')
            
        #%% What part of the EMG do you want to take the percentile of
        # All the EMG data in each succesful trial ('All_EMG')
        # The max EMG in each succesful trial ('Max_EMG')
        norm_method = 'All_EMG'

        #%% Concatenate all the morning & afternoon information
        
        # Time frame
        noon_time_frame = xds_noon.time_frame + xds_morn.time_frame[-1]
        time_frame = np.concatenate((xds_morn.time_frame, noon_time_frame), axis = 0)
        
        # Bin width
        bin_width = xds_morn.bin_width
        
        # EMG
        cat_EMG = np.concatenate((xds_morn.EMG[:,M], xds_noon.EMG[:,M]), axis = 0)
        
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
            
        #%% Pull the EMG corresponding to the extracted time frames
        
        if norm_method == 'All_EMG':
            # Pull all the EMG
            EMG = [[] for ii in range(len(M))]
            for ii in range(len(M)):
                for jj in range(len(total_rewarded_idx)):
                    if jj == 0:
                        EMG[ii] = cat_EMG[int(rewarded_start_idx[jj]) : int(rewarded_end_idx[jj] + 2/bin_width), ii]
                    else:
                        EMG[ii] = np.concatenate((EMG[ii], cat_EMG[int(rewarded_start_idx[jj]) : int(rewarded_end_idx[jj] + 2/bin_width), ii]), axis = 0)
                        
            # Find the percentile of all the EMG
            EMG_Norm_Factor = np.zeros((1,len(M)))
            for ii in range(len(M)):
                EMG_Norm_Factor[0,ii] = np.percentile(EMG[ii], norm_perctile)
                        
        #%% Finding the max EMG per trial
        
        if norm_method == 'Max_EMG':
            max_EMG_pertrial = np.zeros((len(total_rewarded_idx), len(M)))
            for ii in range(len(total_rewarded_idx)):
                for jj in range(len(M)):
                    max_EMG_pertrial[ii,jj] = max(cat_EMG[int(rewarded_start_idx[ii]) : int(rewarded_end_idx[ii] + 2/bin_width), jj])
                    
            # Find the percentile of the max's
            EMG_Norm_Factor = np.zeros((1,len(M)))
            for ii in range(len(M)):
                EMG_Norm_Factor[0,ii] = np.percentile(max_EMG_pertrial[:,ii], norm_perctile)
    
    
        return EMG_Norm_Factor
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    