function con = experimentZero(m)
%experimentZero Default structure for KroneckerBio experimental conditions
%
%   con = experimentZero(m)
% 
%   This function returns a valid, empty experiment. It is used for
%   making experiments constructed by different means consistent so that
%   their structures can be combined.
%
%   Inputs
%   m: [ Kronecker model structure scalar ]
%
%   Outputs:
%   con: [ Kronecker experimental condition structure scalar ]
%       The experiment returned by this function is pretty useless, but the
%       fields can be overwritten. Their meanings are given below.
%       .Type 'Experiment'
%       .Name [ string ]
%           An aribitrary name for the experiment
%       .tF [ nonnegative scalar ]
%           The time at which this experiment ends.
%       .nu [ whole scalar ]
%           The number of inputs
%       .ns [ whole scalar ]
%           The number of seeds
%       .nq [ whole scalar ]
%           The number of input control parameters
%       .nh [ whole scalar ]
%           The number of dose control parameters
%       .s [ nonnegative vector ns ]
%           The values of the seed parameters
%       .q [ real vector nq ]
%           The values of the input control parameters
%       .h [ real vector nh ]
%           The values of the dose control parameters
%       .u [ handle @(t) returns nonnegative matrix nu by numel(t) ]
%           A function handle that returns the value of the inputs at any
%           time t. This function should accept t as a vector.
%       .dudq [ handle @(t) returns real matrix nu by nq ]
%           A function handle that returns the derivative of each input
%           with respect to each input control parameter
%       .d2udq2 [ handle @(t) returns ral matrix nu*nq by nq ]
%           A function handle that returns the derivative of dudq with
%           respect to each input control parameter
%       .d [ handle @(t) returns nonnegative matrix ns by numel(t) ]
%           A function handle that returns the values of the doses at any
%           time t. This function should accept t as a vector.
%       .dddh [ handle @(t) returns real matrix ns by nh ]
%           A function handle that returns the derivative of each dose with
%           respect to each dose control parameter
%       .d2ddh2 [ handle @(t) returns real matrix ns*nh by nh ]
%           A function handle that returns the derivative of dddh with
%           respect to each dose control parameter
%       .inp [ struct scalar ]
%           The input object
%       .dos [ struct scalar ]
%           The dosing object
%       .SteadyState [ logical scalar ]
%           Declares if the system should be run to steady state before
%           the experiment begins
%       .Periodic [ logical scalar ]
%           Declares if the system should be run to a periodic steady state
%           for the duration of the experiment
%       .Discontinuities [ nonegative vector ]
%           Any discontinuous times in the u function must be listed here
%           in order to ensure sucessful integration of the system.
%       .Update [ handle @(s,q,h) returns struct scalar ]
%           This function handle allows the parameter values of the
%           experiment to be changed. Each vector s, q, and h must have the
%           same size as their appropriate counterparts. The structure
%           returned is identical to the one to which Update is attached
%           except that the parameter values have been changed and the
%           appropriate matrices and function handles have been updated.
%       .private [ anything ]
%           Some experimental conditions have extra information associated
%           with them, which is stored here because Matlab does not like
%           stacking structs with different fields.

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Special case of numeric inputs
if isnumeric(m)
    con = emptystruct(m, 'Type', 'Name', 'nu', 'ns', 'nq', 'nh', 's', 'q', 'h', 'u', 'dudq', 'd2udq2', 'd', 'dddh', 'd2ddh2', 'inp', 'dos', 'SteadyState', 'Periodic', 'Discontinuities', 'Update', 'private');
    return
end

con = experimentInitialValue(m, [], [], [], 'Zero');
