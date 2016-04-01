function [con, n_con] = fixCondition(con)
% Standardize experimental conditions structures

assert(is(con, 'Experiment'), 'KroneckerBio:con', 'con must be a vector of experimental conditions')

con = vec(con);
n_con = numel(con);
