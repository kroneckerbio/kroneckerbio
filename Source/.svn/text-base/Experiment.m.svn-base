function con = Experiment(m, tF, x0, steadyState, periodic, u, discontinuities, q, dudq, name)
%Experiment constructs a KroneckerBio experimental conditions structure
%   taking advantage of the full potential of experimental conditions
%
%   con = Experiment(m, tF, x0, steadyState, periodic, u, discontinuities,
%                    q, dudq, name)
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
%   steadyState: [ logical scalar {false} ]
%       Declares if the system should be run to steady state before
%       the experiment begins
%   periodic: [ logical scalar {false} ]
%       Declares if the system should be run to a periodic steady state for
%       the duration of the experiment
%   u: [ handle @(t,q) returns nonnegative matrix nu by numel(t) | 
%       nonnegative vector nu {m.u} ]
%       A function handle that returns the value of the inputs at any time
%       t with associated control parameters q. This function should accept
%       t as a vector.
%   discontinuities: [ nonnegative vector ]
%       Any discontinuous times in the u function must be listed here
%       in order to ensure sucessful integration of the system.
%   q: [ real vector nq {[]} ]
%       The values of the input control parameters
%   dudq: [ handle @(t,q) returns real matrix nu by nq ]
%       A function handle that returns the derivative of each input with
%       respect to each input control parameter. This does not have to be
%       supplied but some Kronecker features will be disabled.
%   name: [ string {''} ]
%       An arbitrary name for the experiment
%
%   Outputs
%   con: [ experiment struct scalar ]
%       The KroneckerBio experimental conditions structure

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 10
    name = [];
    if nargin < 9
        dudq = [];
        if nargin < 8
            q = [];
            if nargin < 7
                discontinuities = [];
                if nargin < 6
                    u = [];
                    if nargin < 5
                        periodic = [];
                        if nargin < 4
                            steadyState = [];
                            if nargin < 3
                                x0 = [];
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Defaults
if isempty(x0)
    x0 = m.x0;
end
if isempty(steadyState)
    steadyState = false;
end
if isempty(periodic)
    periodic = false;
end
if isempty(u) && m.nu > 0
    u = m.u;
end
if isempty(q)
    q = zeros(0,1);
end
if isempty(dudq)
    dudq = zeros(m.nu,0);
end

%% Strip m
temp = m;
clear m;
m.nx = temp.nx;
m.nu = temp.nu;
clear temp

%% Check tF
assert(isscalar(tF) && isreal(tF) && tF >= 0, 'KroneckerBio:Experiment:FinalTime', 'Final time tF must be a scalar real number greater than or equal to zero')

%% Check x0
assert(numel(x0) == m.nx, 'KroneckerBio:Experiment:InitialValueWrongLength', 'Initial values x0 must a length equal to m.nx')

%% Check steadyState
try
    steadyState = logical(steadyState);
catch
    error('KroneckerBio:Experiment:SteadyStateType', 'steadyState must be true or false, or convertible to such')
end

%% Check periodic
try
    periodic = logical(periodic);
catch
    error('KroneckerBio:Experiment:PeriodicType', 'periodic must be true or false, or convertible to such')
end

%% Standardize q as a numeric column vector
assert(isnumeric(q), 'KroneckerBio:Experiment:ControlType', 'Control parameters q must be a vector of numerics')
q = vec(q);

%% Standardize dudq is possible
if isempty(q)
    dudq = @(t,q)zeros(m.nu,0);
end
if isnumeric(u)
    dudq = @(t,q)zeros(m.nu,numel(q));
end

%% Standardize different types of inputs as @(t,q)
% Numeric
if isnumeric(u)
    value = vec(u);
    assert(numel(u) == m.nu, 'KroneckerBio:ExperimentBasic:InputWrongLength', 'Numeric inputs u must have a length of m.nu')
    u = @(t,q)piecewisestep(t, [-1; 0], [zeros(1,numel(value)); value.']).';
    dudq = @(t,q)zeros(m.nu, numel(q));
end

% Function handle as @(t)
if isa(u, 'function_handle') && nargin(u) == 1
    u = @(t,q)u(t);
end

assert(isa(u, 'function_handle') && nargin(u) == 2, 'KroneckerBio:Experiment:InputType', 'Input (u) must be a positive vector nu or a function handle @(t) or @(t,q) returning a positive vector nu')

%% Standardize discontinuities vector as a unique, sorted, column vector
assert(isnumeric(discontinuities) && all(discontinuities >= 0), 'KroneckerBio:Experiment:DiscontinuitiesType', 'discontinuities must be a positive vector')
% Empty
if isempty(discontinuities)
    discontinuities = zeros(0,1);
end

discontinuities = unique(discontinuities);

%% Standardize name as a string
if isempty(name)
    name = '';
end
assert(ischar(name), 'KroneckerBio:Experiment:NameType', 'name must be a string')

%% Build experiment
con.Type = 'Experiment';
con.Name = name;
con.tF = tF;
con.x0 = x0;
con.u  = @(t)u(t,q);
con.q  = q;
con.dudq = @(t)dudq(t,q);
con.nq = numel(q);
con.SteadyState = steadyState;
con.Periodic = periodic;
con.Discontinuities = discontinuities;
con.Update = @Update;

    function conOut = Update(x0, q)
        conOut = Experiment(m, tF, x0, steadyState, periodic, u, discontinuities, q, dudq, name);
    end
end