# -*- coding: utf-8 -*-
from EMGIndex import EMG_Index
import numpy as np

class SingleSession_EMGZero():
    
    def Single_Session_EMG_Zero(xds, muscle_groups, zero_method, zero_EMG):
        
        #%% Find the EMG index
        
        M = EMG_Index.EMG_Indexing(xds, muscle_groups)
        
        #%% End the function if you are not zeroing the EMG
        if zero_EMG == 0:
            print('EMG Will Not Be Zeroed')
            EMG_Zero_Factor = np.zeros([1,len(M)])
            return EMG_Zero_Factor
        else:
            print('EMG Will Be Zeroed')
            
        #%% Extract the EMG
        EMG = np.zeros([len(xds.time_frame), len(M)])
        for ii in range(len(M)):
            EMG[:,ii] = xds.EMG[:,M[ii]]
            
        #%% Run the moving average
        if zero_method == 'Window':
            
            # Moving average window size (200 ms)
            window_size = 200
            
            # Calculate the moving average
            move_avg = np.zeros([len(xds.time_frame) - window_size + 1, len(M)])
            for ii in range(len(M)):
                for jj in range(len(EMG) - window_size + 1):
                    move_avg[jj,ii] = np.mean(EMG[jj:jj+window_size, ii])
                    
            # Find the minimum for each muscle
            EMG_Zero_Factor = np.zeros([1,len(M)])
            for ii in range(len(M)):
                EMG_Zero_Factor[0,ii] = min(move_avg[:,ii])
                
        #%% Find the 5th percentile
        if zero_method == 'Percentile':
            
            EMG_Zero_Factor = np.zeros([1,len(M)])
            for ii in range(len(M)):
                EMG_Zero_Factor[0,ii] = np.percentile(EMG[:,ii], 5)
        
                
        return EMG_Zero_Factor            
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    