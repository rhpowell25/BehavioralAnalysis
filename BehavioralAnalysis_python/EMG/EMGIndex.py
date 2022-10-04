# -*- coding: utf-8 -*-

import numpy as np

class EMG_Index:
    
    def EMG_Indexing(xds, muscle_groups):
        
        #%% Define the muscles
        
        muscle_names = []
        
        if muscle_groups == 'Flex':
            muscle_names.append('FCR')
            muscle_names.append('FCU')
            
        if muscle_groups == 'Exten':
            muscle_names.append('ECR')
            muscle_names.append('ECU')
            
        if muscle_groups == 'Both':
            muscle_names.append('FCR')
            muscle_names.append('FCU')
            muscle_names.append('ECR')
            muscle_names.append('ECU')
            
        if muscle_groups == 'Uln_Dev':
            muscle_names.append('FCU')
            muscle_names.append('ECU')
            
        if muscle_groups == 'Rad_Dev':
            muscle_names.append('FCR')
            muscle_names.append('ECR')
        
        if muscle_groups == 'Grasp':
            muscle_names.append('FDP')
            muscle_names.append('FDS')
        
        if muscle_groups == 'Custom':
            muscle_names.append('FDS1')
        
        # If muscle_names is still blank
        if muscle_names == []:
            muscle_names.append(muscle_groups)
            
        # Find the indices of the muscles of interest
        muscle_idx = []
        cc = 0
        for ii in range(len(muscle_names)):
            muscle_idx.append([idx for idx, s in enumerate(xds.EMG_names) if muscle_names[ii] in s])
            # Find how many muscles there are
            cc = cc + len(muscle_idx[ii])
            
        # Concatenate the indices
        M = []
        for ii in range(len(muscle_idx)):
            for jj in range(len(muscle_idx[ii])):
                M.append(muscle_idx[ii][jj])
            
        if muscle_groups == 'All':
            M = list(range(len(xds.EMG_names)))
            
            
        return M
    
    
    