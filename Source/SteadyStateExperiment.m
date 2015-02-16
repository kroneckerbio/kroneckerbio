function con = SteadyStateExperiment(m, tF, s, u, d, discontinuities, q, dudq, dddq, name)
%Experiment constructs a KroneckerBio experimental conditions structure
%   taking advantage of the full potential of experimental conditions
%
%   con = Experiment(m, tF, s, steadyState, periodic, u, d,
%                    discontinuities, q, dudq, dddq, name)
%
%   The inputs to this function allow one to set all the variables that are
%   permitted on a KroneckerBio experimental conditions structure.
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model these experiments will be run for
%   tF: [ nonnegative scalar ]
%       The time at which this experiment ends.
%   s: [ nonnegative vector ns {m.s} ]
%       The values of the seed parameters
%   u: [ handle @(t,q) returns nonnegative matrix nu by numel(t) | 
%       nonnegative vector nu {m.u} ]
%       A function handle that returns the value of the inputs at any time
%       t with associated control parameters q. This function should accept
%       t as a vector.
%   d: [ handle @(t,q) returns nonnegative matrix ns by numel(t) |
%        nonegative vector ns {@(t,q)zeros(m.ns,0)} ]
%       A function that returns the values of the doses at any time t
%       sharing control parameters with u. Doses are applied via the seeds
%       to instantaneously update the state of the system. Must defined for
%       all t and return zeros(ns,1) for times at which a dose is not
%       given. If d is a numeric vector, then that dose is applied at all
%       times listed in discontinuities.
%   discontinuities: [ nonnegative vector {[]} ]
%       Any dosing times in the d function and any discontinuous times in
%       the u function must be listed here in order to ensure successful
%       evaluation of d and u.
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
%
%   For the meanings of the fields of con see "help Uzero"

% (c) 2014 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 12
    name = [];
    if nargin < 11
        dddq = [];
        if nargin < 10
            dudq = [];
            if nargin < 9
                q = [];
                if nargin < 8
                    discontinuities = [];
                    if nargin < 7
                        d = [];
                        if nargin < 6
                            u = [];
                            if nargin < 5
                                periodic = [];
                                if nargin < 4
                                    steadyState = [];
                                    if nargin < 3
                                        s = [];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Check m
assert(isscalar(m), 'KroneckerBio:Experiment:ScalarModel', 'Model m must be scalar')

%% Defaults
if isempty(s)
    s = m.s;
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
if isempty(d)
    d = @(t,q)zeros(m.ns,1);
end
if isempty(q)
    q = zeros(0,1);
end

%% Strip m
temp = m;
clear m;
m.s = temp.s;
m.ns = temp.ns;
m.nu = temp.nu;
clear temp

%% Check tF
assert(isscalar(tF) && isreal(tF) && tF >= 0, 'KroneckerBio:Experiment:FinalTime', 'Final time tF must be a scalar real number greater than or equal to zero')

%% Check x0
assert(numel(s) == m.ns, 'KroneckerBio:Experiment:InitialValueWrongLength', 'Seed parameters s must a length equal to m.ns')

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

if periodic; warning('KroneckerBio:Experiment:PeriodicNotImplemented', 'The Periodic option has not yet been implemented in KroneckerBio and will be assumed false'); periodic = false; end

%% Standardize q as a numeric column vector
assert(isnumeric(q), 'KroneckerBio:Experiment:ControlType', 'Control parameters q must be a vector of numerics')
q = vec(q);

%% Standardize discontinuities vector as a unique, sorted, column vector
assert(isnumeric(discontinuities) && all(discontinuities >= 0), 'KroneckerBio:Experiment:DiscontinuitiesType', 'discontinuities must be a nonnegative vector')
% Empty
if isempty(discontinuities)
    discontinuities = zeros(0,1);
end

discontinuities = vec(unique(discontinuities));

%% Infer dudq and dddq if possible
if isempty(q)
    dudq = @(t,q)zeros(m.nu,0);
    dddq = @(t,q)zeros(m.ns,0);
end

%% Standardize different types of inputs as @(t,q)
% Numeric
if isnumeric(u)
    value_u = vec(u);
    assert(numel(u) == m.nu, 'KroneckerBio:ExperimentBasic:InputWrongLength', 'Numeric inputs u must have a length of m.nu')
    u = @(t,q)piecewisestep(t, [-1; 0], [zeros(1,numel(value_u)); value_u.']).';
    dudq = @(t,q)zeros(m.nu, numel(q));
end

% Function handle as @(t)
if isa(u, 'function_handle') && nargin(u) == 1
    u = @(t,q)u(t);
    dudq = @(t,q)zeros(m.nu, numel(q));
end

assert(isa(u, 'function_handle') && nargin(u) == 2, 'KroneckerBio:Experiment:InputType', 'Input (u) must be a positive vector nu or a function handle @(t) or @(t,q) returning a positive vector nu')

%% Standardize different types of doses as @(t,q)
% Numeric
if isnumeric(d)
    assert(numel(d) == m.ns, 'KroneckerBio:ExperimentBasic:DoseWrongLength', 'Numeric doses d must have a length of m.ns')
    d = makeDiscreteDoses(vec(d), discontinuities);
    dddq = @(t,q)zeros(m.nu, numel(q));
end

% Function handle as @(t)
if isa(d, 'function_handle') && nargin(d) == 1
    d = @(t,q)d(t);
    dddq = @(t,q)zeros(m.nu, numel(q));
end

assert(isa(d, 'function_handle') && nargin(d) == 2, 'KroneckerBio:Experiment:DosingType', 'Dosing (d) must be a positive vector ns or a function handle @(t) or @(t,q) returning a positive vector ns')

%% Standardize name as a string
if isempty(name)
    name = '';
end
assert(ischar(name), 'KroneckerBio:Experiment:NameType', 'name must be a string')

%% Build experiment
con.Type = 'Experiment';
con.Name = name;
con.tF = tF;
con.s = s;
con.u  = @(t)u(t,q);
con.d  = @(t)d(t,q);
con.q  = q;
con.dudq = @(t)dudq(t,q);
con.dddq = @(t)dddq(t,q);
con.nq = numel(q);
con.SteadyState = steadyState;
con.Periodic = periodic;
con.Discontinuities = discontinuities;
con.Update = @update;

    function conOut = update(s, q)
        conOut = Experiment(m, tF, s, steadyState, periodic, u, d, discontinuities, q, dudq, dddq, name);
    end
end

function handle = makeDiscreteDoses(d, discontinuities)

handle = @discreteDoses;

    function s = discreteDoses(t,q)
        if any(t == discontinuities)
            s = d;
        else
            s = zeros(numel(d),1);
        end
    end
end