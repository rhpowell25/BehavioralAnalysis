function [Sigma_Force] = Sum_Force(Task, Force)

%% Sum the two force transducers

Sigma_Force = struct([]);
for ii = 1:length(Force)
    if strcmp(Task, 'WB')
        % Vector sum
        Sigma_Force{ii,1} = sqrt(Force{ii,1}(:,1).^2 + Force{ii,1}(:,2).^2);
    else
        % Algebraic sum
        Sigma_Force{ii,1} = Force{ii,1}(:,1) + Force{ii,1}(:,2);
    end
end








