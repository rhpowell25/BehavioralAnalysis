# -*- coding: utf-8 -*-

from EMGIndex import EMG_Index
import numpy as np

class TrialLength_RxnTime:    
    
    def __init__(self, xds):
        
        #%% Basic settings, some variable extractions, & definitions
        
        # Define the period before Go-Cue & end to measure
        time_before_gocue = 0.4
        
        # Find the EMG index
        if xds._lab_data__meta['task_name'] == 'PG':
            muscle_groups = 'Grasp'
            M = EMG_Index.EMG_Indexing(xds, muscle_groups)
        
        #%% Indexes for rewarded trials in all directions
        
        # Extract the trial directions
        target_dir_idx = np.round(xds.trial_target_dir)
        
        # Select the first direction (start with the minimum direction value)
        target_dirs = np.unique(target_dir_idx)
        # Count the number of directions used
        num_dir = len(target_dirs)
        
        # Set the loop counter
        cc = 0
        
        #%% Begin the loop through all directions
        for jj in range(num_dir):
            
            #%% Find the EMG index
            if xds._lab_data__meta['task_name'] == 'WS':
                if target_dirs[jj] == 0 and xds._lab_data__meta['hand'] == 'Left':
                    muscle_groups = 'Flex'
                if target_dirs[jj] == 90:
                    muscle_groups = 'Rad_Dev'
                if target_dirs[jj] == 180 and xds._lab_data__meta['hand'] == 'Left':
                    muscle_groups = 'Exten'
                if target_dirs[jj] == -90:
                    muscle_groups = 'Uln_Dev'   
                M = EMG_Index.EMG_Indexing(xds, muscle_groups)
                
            #%% Indexes for rewarded trials
            total_rewarded_idx = np.argwhere(np.in1d(xds.trial_result == 'R', xds.trial_target_dir == target_dirs[jj])).reshape(-1,)
            
            #%% Find the number of targets in that particular direction
            try:
                tgt_Center_idx = xds.trial_info_table_header.index('tgtCenter')
            except ValueError:
                tgt_Center_idx = xds.trial_info_table_header.index('tgtCtr')
            
            # Pull the target center coordinates of each succesful trial
            tgt_cntrs = []
            for ii in range(len(total_rewarded_idx)):
                tgt_cntrs.append(xds.trial_info_table[tgt_Center_idx][total_rewarded_idx[ii]])
            
            # Convert the the cartesian coordinates into polar coordinates
            target_cntrs = np.zeros(len(total_rewarded_idx))
            for ii in range(len(total_rewarded_idx)):
                target_cntrs[ii] = (tgt_cntrs[ii][0][0]**2 + tgt_cntrs[ii][1][0]**2)**0.5
            
            # Confirm both sessions use consistent target centers
            unique_targets = np.unique(target_cntrs)
            
            if jj == 0:
                self.target_centers = np.zeros(num_dir*len(unique_targets))
                self.target_dirs = np.zeros(num_dir*len(unique_targets))
                self.rxn_time = [[] for ii in range(len(unique_targets))]
                self.trial_length = [[] for ii in range(len(unique_targets))]
                
            # Remove output variables when there are no unique targets
            if len(unique_targets) == 0:
                self.target_centers[cc] = []
                self.target_dirs[cc] = []
            
            #%% Redefine the rewarded_idx according to the target center
            for kk in range(len(unique_targets)):
                
                #%% Indexes for rewarded trials
                rewarded_idx = total_rewarded_idx[np.asarray(np.where(target_cntrs == unique_targets[kk])).reshape(-1,)]
                
                #%% Loops to extract only rewarded trials
                # End times for succesful trials
                rewarded_end_time = np.zeros(len(rewarded_idx))
                for ii in range(len(rewarded_idx)):
                    rewarded_end_time[ii] = xds.trial_end_time[rewarded_idx[ii]]
                    
                # Go-cue's for succesful trials
                rewarded_gocue_time = np.zeros(len(rewarded_idx))
                for ii in range(len(rewarded_idx)):
                    rewarded_gocue_time[ii] = xds.trial_gocue_time[rewarded_idx[ii]]
                    
                #%% Skip the function if less than 5 succesful trials in that direction
        
                if len(rewarded_idx) < 5:
                    self.target_dirs[cc] = [];
                    self.target_centers[cc] = [];
                    continue
            
                #%% Extracting EMG & time during successful trials
                
                # EMG & time measured during each succesful trial
                Baseline_EMG = [[] for ii in range(len(rewarded_idx))] # EMG during each succesful trial
                Baseline_timings = [] # Time points during each succesful trial
                for ii in range(len(rewarded_idx)):
                    Baseline_idx = np.where((xds.time_frame > rewarded_gocue_time[ii] - time_before_gocue)*
                         (xds.time_frame < rewarded_gocue_time[ii]))
                    for mm in range(len(M)):
                        if mm == 0:
                            Baseline_EMG[ii] = xds.EMG[Baseline_idx, M[mm]]
                        else:
                            Baseline_EMG[ii] = np.vstack((Baseline_EMG[ii], xds.EMG[Baseline_idx, M[mm]]))
                    Baseline_EMG[ii] = np.transpose(Baseline_EMG[ii])
                    Baseline_timings.append(xds.time_frame[Baseline_idx])
                    
                # EMG & time measured during each succesful trial
                Trial_EMG = [[] for ii in range(len(rewarded_idx))] # EMG during each suxxesful trial
                Trial_timings = [] # Time points during each succesful trial
                for ii in range(len(rewarded_idx)):
                    Trial_idx = np.where((xds.time_frame > rewarded_gocue_time[ii])*
                         (xds.time_frame < rewarded_end_time[ii]))
                    for mm in range(len(M)):
                        if mm == 0:
                            Trial_EMG[ii] = xds.EMG[Trial_idx, M[mm]]
                        else:
                            Trial_EMG[ii] = np.vstack((Trial_EMG[ii], xds.EMG[Trial_idx, M[mm]]))
                    Trial_EMG[ii] = np.transpose(Trial_EMG[ii])
                    Trial_timings.append(xds.time_frame[Trial_idx])
                    
                #%% Find the mean of the EMG
                mean_Baseline_EMG = [[] for ii in range(len(Baseline_EMG))]
                for ii in range(len(Baseline_EMG)):
                    for mm in range(len(Baseline_EMG[ii])):
                        mean_Baseline_EMG[ii].append(np.mean(Baseline_EMG[ii][mm]))
            
                mean_Trial_EMG = [[] for ii in range(len(Trial_EMG))]
                for ii in range(len(Trial_EMG)):
                    for mm in range(len(Trial_EMG[ii])):
                        mean_Trial_EMG[ii].append(np.mean(Trial_EMG[ii][mm]))
                
                #%% Find the standard deviations of the baseline EMG
                    std_Baseline_EMG = np.zeros(len(mean_Baseline_EMG))
                    for ii in range(len(mean_Baseline_EMG)):
                        std_Baseline_EMG[ii] = np.std(mean_Baseline_EMG[ii])
            
                #%% Find the trial EMG that exceeds 2 std of the baseline EMG
                    rxn_time_EMG = np.zeros(len(mean_Trial_EMG))
                    rxn_time_idx = np.zeros(len(mean_Trial_EMG))
                    for ii in range(len(mean_Trial_EMG)):
                        rxn_time_EMG_idx = np.where(mean_Trial_EMG[ii] >= 
                            np.mean(mean_Baseline_EMG[ii]) + 2*std_Baseline_EMG[ii])[0]
                        if len(rxn_time_EMG_idx) == 0: # If nothing exceeds 2 std try 1 std
                            rxn_time_EMG_idx = np.where(mean_Trial_EMG[ii] >= 
                                np.mean(mean_Baseline_EMG[ii]) + std_Baseline_EMG[ii])[0]
                            if len(rxn_time_EMG_idx) == 0: # if nothing exceeds 1 std just find the max of that trial
                                rxn_time_EMG_idx = np.where(mean_Baseline_EMG[ii] ==
                                    np.max(mean_Baseline_EMG[ii]))[0]
                        rxn_time_idx[ii] = rxn_time_EMG_idx[0]
                        rxn_time_EMG[ii] = Trial_timings[ii][int(rxn_time_idx[ii])]
                                
                #%% Defining the output variables
                # Target Direction
                self.target_dirs[cc] = target_dirs[jj]
                # Target Centers
                self.target_centers[cc] = round(unique_targets[kk])
            
                # Reaction time
                self.rxn_time[cc] = rxn_time_EMG - rewarded_gocue_time
                # Trial Length
                self.trial_length[cc] = rewarded_end_time - rewarded_gocue_time
            
                # Adding to the loop counter
                cc = cc + 1

















    
    
    
    