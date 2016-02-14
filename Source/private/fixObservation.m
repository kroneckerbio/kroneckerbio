function [obs, n_obs] = fixObservation(obs, n_con)
% Standardize observation scheme structures

if isnumeric(obs)
    obs = observationAll(obs);
end

%assert(is(obs, 'Observation'), 'KroneckerBio:obs', 'obs must be a matrix of observation schemes')
assert(size(obs,2) == n_con, 'KroneckerBio:obsSize', 'obs must have a number of columns equal to numel(con)')

n_obs = size(obs,1);
