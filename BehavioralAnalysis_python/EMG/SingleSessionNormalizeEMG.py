# -*- coding: utf-8 -*-
from EMGIndex import EMG_Index
import numpy as np

class SingleSession_NormalizeEMG():
    
    def Single_Session_Normalize_EMG(xds, muscle_groups, norm_perctile, norm_EMG):
        
        #%% Find the EMG index
        
        M = EMG_Index.EMG_Indexing(xds, muscle_groups)
        
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
        
        #%% If normalizing to the entire experiment
        if xds.has_trials == 0:
            
            #%% Find the percentile of all the EMG
            
            if norm_method == 'All_EMG':
                
                EMG_Norm_Factor = np.zeros((1, len(M)))
                for ii in range(len(M)):
                    EMG_Norm_Factor[0,ii] = np.percentile(xds.EMG[:,M[ii]], norm_perctile)
                    
            #%% Finding the max EMG per trial
            
            if norm_method == 'Max_EMG':
                
                max_EMG_perchannel = np.zeros((1, len(M)))
                for jj in range(len(M)):
                    max_EMG_perchannel[jj] = max(xds.EMG[:,M[jj]])
                    
                # Find the percentile of the max's
                EMG_Norm_Factor = np.percentile(max_EMG_perchannel, norm_perctile)

        #%% If normalizing to only succesful trials
        if xds.has_trials == 1:
            
            # Bin width
            bin_width = xds.bin_width
            
            #%% Index for rewarded trials
            
            total_rewarded_idx = np.argwhere(xds.trial_results == 'R').reshape(-1,)
            
            #%% Loops to extract only rewarded trials
            
            # Rewarded start times
            rewarded_start_time = np.zeros(len(total_rewarded_idx))
            for ii in range(len(total_rewarded_idx)):
                rewarded_start_time[ii] = xds.trial_start_time[total_rewarded_idx[ii]]
                
            # Rewarded go-cue's
            rewarded_gocue_time = np.zeros(len(total_rewarded_idx))
            for ii in range(len(total_rewarded_idx)):
                rewarded_gocue_time[ii] = xds.trial_gocue_time[total_rewarded_idx[ii]]
                    
            # Rewarded end times
            rewarded_end_time = np.zeros(len(total_rewarded_idx))
            for ii in range(len(total_rewarded_idx)):
                rewarded_end_time[ii] = xds.trial_end_time[total_rewarded_idx[ii]]
                
            #%% Pulling the timeframe of the succesful trials
            # Find the rewarded start times in the whole trial time frame
            rewarded_start_idx = np.zeros(len(total_rewarded_idx))
            for ii in range(len(total_rewarded_idx)):
                rewarded_start_idx[ii] = np.asarray(np.where(xds.time_frame == rewarded_start_time[ii]))
                
            # Find the rewarded end times in the whole trial time frame
            rewarded_end_idx = np.zeros(len(total_rewarded_idx))
            for ii in range(len(total_rewarded_idx)):
                rewarded_end_idx[ii] = np.asarray(np.where(xds.time_frame == rewarded_end_time[ii]))  
                
            #%% Pull the EMG corresponding to the extracted time frames
            
            if norm_method == 'All_EMG':
                # Pull all the EMG
                EMG = [[] for ii in range(len(M))]
                for ii in range(len(M)):
                    for jj in range(len(total_rewarded_idx)):
                        if jj == 0:
                            EMG[ii] = xds.EMG[int(rewarded_start_idx[jj]) : int(rewarded_end_idx[jj] + 2/bin_width), M[ii]]
                        else:
                            EMG[ii] = np.concatenate((EMG[ii], xds.EMG[int(rewarded_start_idx[jj]) : int(rewarded_end_idx[jj] + 2/bin_width), M[ii]]), axis = 0)
                            
                # Find the percentile of all the EMG
                EMG_Norm_Factor = np.zeros((1,len(M)))
                for ii in range(len(M)):
                    EMG_Norm_Factor[0,ii] = np.percentile(EMG[ii], norm_perctile)
                            
            #%% Finding the max EMG per trial
            
            if norm_method == 'Max_EMG':
                max_EMG_pertrial = np.zeros((len(total_rewarded_idx), len(M)))
                for ii in range(len(total_rewarded_idx)):
                    for jj in range(len(M)):
                        max_EMG_pertrial[ii,jj] = max(xds.EMG[int(rewarded_start_idx[ii]) : int(rewarded_end_idx[ii] + 2/bin_width), M[jj]])
                        
                # Find the percentile of the max's
                EMG_Norm_Factor = np.zeros((1,len(M)))
                for ii in range(len(M)):
                    EMG_Norm_Factor[0,ii] = np.percentile(max_EMG_pertrial[:,ii], norm_perctile)
    
    
        return EMG_Norm_Factor
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    