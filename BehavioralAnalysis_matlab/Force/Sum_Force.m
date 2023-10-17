function [Sigma_Force] = Sum_Force(Task, Force)

%% Sum the two force transducers

if iscell(Force)
    Sigma_Force = struct([]);
    for ii = 1:length(Force)
        if strcmp(Task, 'WB')
            % Vector sum
            Sigma_Force{ii,1} = sqrt(Force{ii,1}(:,1).^2 + Force{ii,1}(:,2).^2);
        else
            % Algebraic sum
            Sigma_Force{ii,1} = sum(Force{ii,1}, 2);
        end
    end
else
    if strcmp(Task, 'WB')
        % Vector sum
        Sigma_Force = sqrt(Force(:,1).^2 + Force(:,2).^2);
    else
        % Algebraic sum
        Sigma_Force = sum(Force, 2);
    end
end








