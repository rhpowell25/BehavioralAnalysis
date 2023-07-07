function [muscle_groups] = Hand_Muscle_Match(xds, target_dir)

%% Define the muscle groups interest, based on the task & the unit's preferred direction

if contains(xds.meta.rawFileName, 'PG')
    muscle_groups = 'Grasp';
end

if isequal(xds.meta.task, 'WS') || isequal(xds.meta.task, 'WB')
    if isequal(target_dir, 0)
        if strcmp(xds.meta.hand, 'Left')
            muscle_groups = 'Flex';
        end
        if strcmp(xds.meta.hand, 'Right')
            muscle_groups = 'Ext';
        end
    end
    if isequal(target_dir, 45)
        if strcmp(xds.meta.hand, 'Left')
            muscle_groups = 'FCR';
        end
        if strcmp(xds.meta.hand, 'Right')
            muscle_groups = 'ECR';
        end
    end
    if isequal(target_dir, 90)
        muscle_groups = 'Rad_Dev';
    end
    if isequal(target_dir, 135) 
        if strcmp(xds.meta.hand, 'Left')
            muscle_groups = 'ECR';
        end
        if strcmp(xds.meta.hand, 'Right')
            muscle_groups = 'FCR';
        end
    end
    if isequal(target_dir, 180) 
        if strcmp(xds.meta.hand, 'Left')
            muscle_groups = 'Ext';
        end
        if strcmp(xds.meta.hand, 'Right')
            muscle_groups = 'Flex';
        end
    end
    if isequal(target_dir, -135) 
        if strcmp(xds.meta.hand, 'Left')
            muscle_groups = 'ECU';
        end
        if strcmp(xds.meta.hand, 'Right')
            muscle_groups = 'FCU';
        end
    end
    if isequal(target_dir, -90)
        muscle_groups = 'Uln_Dev';
    end
    if isequal(target_dir, -45) 
        if strcmp(xds.meta.hand, 'Left')
            muscle_groups = 'FCU';
        end
        if strcmp(xds.meta.hand, 'Right')
            muscle_groups = 'ECU';
        end
    end
end





