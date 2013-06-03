function con = ExperimentBasic(m, tF, x0, u, discontinuities, name)
%ExperimentBasic constructs a KroneckerBio experimental conditions
%   structure taking advantage of the full potential of experimental
%   conditions
%
%   con = ExperimentBasic(m, tF, x0, u, discontinuities, name)
%
%   The inputs to this function allow one to set all the variables that are
%   permitted on a KroneckerBio experimental conditions structure.
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model these experiments will be run for
%   tF: [ nonnegative scalar ]
%       The time at which this experiment ends.
%   x0: [ nonnegative vector nx {m.x0} ]
%       The values of the initial conditions of each state species
%   u: [ handle @(t) returns nonnegative matrix nu by numel(t) | 
%       nonnegative vector nu {m.u} ]
%       A function handle that returns the value of the inputs at any time
%       t. This function should accept t as a vector.
%   discontinuities: [ nonnegative vector ]
%       Any discontinuous times in the u function must be listed here
%       in order to ensure sucessful integration of the system.
%   name: [ string {''} ]
%       An arbitrary name for the experiment
%
%   Outputs
%   con: [ experiment struct scalar ]
%       The KroneckerBio experimental conditions structure

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Sanitize inputs
if nargin < 6
    name = 'UnnamedExperiment';
    if nargin < 5
        discontinuities = [];
    end
end

%% Strip m
temp = m;
clear m;
m.nx = temp.nx;
m.nu = temp.nu;
m.nq = temp.nq;
clear temp

%% Check tF
assert(isscalar(tF) && isreal(tF) && tF >= 0, 'KroneckerBio:ExperimentBasic:FinalTime', 'Final time tF must be a scalar real number greater than or equal to zero')

%% Check x0
assert(numel(x0) == m.nx, 'KroneckerBio:ExperimentBasic:InitialValueWrongLength', 'Initial values x0 must a length equal to m.nx')

%% Standardize different types of inputs as @(t)
% Empty
if isempty(u)
    assert(m.nu == 0, 'KroneckerBio:ExperimentBasic:EmptyInput', 'Input must be specified for ExperimentBasic')
end

% Numeric
if isnumeric(u)
    value = vec(u);
    % Standardize empty as zeros(0,1)
    if isempty(value)
        value = zeros(0,1);
    end
    assert(numel(u) == m.nu, 'KroneckerBio:ExperimentBasic:InputWrongLength', 'Numeric inputs u must have a length of m.nu')
    u = vec(u);
    u = @(t)repmat(value, 1,numel(t));
end

% Function handle as @(t,q)
if isa(u, 'function_handle') && nargin(u) == 2
    u = @(t)u(t,zeros(0,1));
end

%% Standardize discontinuities vector as a unique, sorted, column vector
% Empty
if isempty(discontinuities)
    discontinuities = zeros(0,1);
end

discontinuities = unique(discontinuities);

%% Standardize name as a string
if isempty(name)
    name = '';
end

%% Build experiment
con.Type = 'Experiment';
con.Name = name;
con.tF = tF;
con.x0 = x0;
con.u  = u;
con.q  = zeros(0,1);
con.dudq = @(t)zeros(m.nu,0);
con.nq = 0;
con.SteadyState = false;
con.Periodic = false;
con.Discontinuities = discontinuities;
con.Update = @Update;

    function conOut = Update(x0, q)
        conOut = ExperimentBasic(m, tF, x0, u, discontinuities, name);
    end
end