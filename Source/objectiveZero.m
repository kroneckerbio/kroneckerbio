function obj = Gzero(dims)
%Gzero Default structure for Kronecker Bio objective functions
% 
%   obj = Gzero(size)
%
%   This function is a blank objective function. It is used to fill in the
%   fields with appropriately zero defaults when they are not provided.
%
%   Inputs
%   dims: [ nonnegative integer vector ]
%       The size of the blank objective structure array
%
%   Outputs:
%   obj: [ Kronecker objective structure scalar ]
%       The meaning of the fields in a Kronecker objective structure is
%       given below.
%       .Type [ 'Objective' 'Objective.Information' 'Objective.Data' ]
%           The Type field can be used to see what analyses can be done on
%           the objective structure.
%           All objective structures must allow for fitting an objective
%           function. Objective structures with only this ability will be
%           labeled with 'Objective'. For some objective structures, it is
%           possible to apply information theory to the objective
%           structure. These objectives will be labeled with
%           'Objective.Information'. If, in addition to information theory,
%           the objective is derived from data, the objective will be
%           labeled 'Objective.Data'.
%       .Name [ string ]
%           An aribitrary name for the objective
%       .Continuous [ logical scalar ]
%           Declares if this objective has a continuous objective function
%       .Complex [ logical scalar ]
%           Declares if this objective requires the entire simulation in
%           order to evaluate the objective function
%       .DiscreteTimes [ nonnegative vector ]
%           If Complex is false, then all discrete time points at which the
%           objective function must be evaluated must be put here.
%       .g [ handle @(t,x,u) returns real scalar ]
%           This is the continuous objective function. This function
%           returns the derivative of the goal function with respect to
%           time. The goal value comes from the numerical integration of
%           this function.
%       .dgdx [ handle @(t,x,u) returns real vector nx ]
%           The partial derivative of g wrt x
%       .d2gdx2 [ handle @(t,x,u) returns real matrix nx by nx ]
%           The partial derivative of dgdx wrt x
%       .G [ handle @(sol) returns real scalar ]
%           This is the discrete objective function. This function
%           evaluates the simulation of the system and returns the
%           objective value. If Complex is true for any objective structure
%           attached to an experiment, sol will be an ODE solution
%           evaluatable by deval(). Otherwise, sol.x will include
%           DiscreteTimes and the corresponding simulation values will be
%           in sol.y.
%       .dGdx [ handle @(t,sol) returns real vector nx ]
%           The derivative of G wrt x at the scalar time point t.
%       .dGdk [ handle @(t,sol) returns real vector nk ]
%           The derivative of G wrt k at the scalar time point t.
%       .dGds [ handle @(t,sol) returns real vector ns ]
%           The derivative of G wrt s at the scalar time point t.
%       .dGdq [ handle @(t,sol) returns real vector nq ]
%           The derivative of G wrt q at the scalar time point t.
%       .dGdh [ handle @(t,sol) returns real vector nh ]
%           The derivative of G wrt h at the scalar time point t.
%       .d2Gdx2 [ handle @(t,sol) returns real vector nx by nx ]
%           The derivative of dGdx wrt x at the scalar time point t.
%       .d2Gdk2 [ handle @(t,sol) returns real vector nk by nk ]
%           The derivative of dGdk wrt k at the scalar time point t.
%       .d2Gds2 [ handle @(t,sol) returns real vector ns by ns ]
%           The derivative of dGds wrt s at the scalar time point t.
%       .d2Gdq2 [ handle @(t,sol) returns real vector nq by nq ]
%           The derivative of dGdq wrt q at the scalar time point t.
%       .d2Gdh2 [ handle @(t,sol) returns real vector nh by nh ]
%           The derivative of dGdq wrt h at the scalar time point t.
%       .d2Gdkdx [ handle @(t,sol) returns real vector nx by nk ]
%           The derivative of dGdx wrt k at the scalar time point t.
%       .d2Gdsdx [ handle @(t,sol) returns real vector nx by ns ]
%           The derivative of dGdx wrt s at the scalar time point t.
%       .d2Gdqdx [ handle @(t,sol) returns real vector nx by nq ]
%           The derivative of dGdx wrt q at the scalar time point t.
%       .d2Gdhdx [ handle @(t,sol) returns real vector nx by nh ]
%           The derivative of dGdx wrt h at the scalar time point t.
%       .d2Gdxdk [ handle @(t,sol) returns real vector nk by nx ]
%           The derivative of dGdk wrt x at the scalar time point t.
%       .d2Gdsdk [ handle @(t,sol) returns real vector nk by s ]
%           The derivative of dGdk wrt s at the scalar time point t.
%       .d2Gdqdk [ handle @(t,sol) returns real vector nk by nq ]
%           The derivative of dGdk wrt q at the scalar time point t.
%       .d2Gdhdk [ handle @(t,sol) returns real vector nk by nh ]
%           The derivative of dGdk wrt h at the scalar time point t.
%       .d2Gdxds [ handle @(t,sol) returns real vector ns by nx ]
%           The derivative of dGds wrt x at the scalar time point t.
%       .d2Gdkds [ handle @(t,sol) returns real vector ns by nk ]
%           The derivative of dGds wrt k at the scalar time point t.
%       .d2Gdqds [ handle @(t,sol) returns real vector ns by ns ]
%           The derivative of dGds wrt q at the scalar time point t.
%       .d2Gdhds [ handle @(t,sol) returns real vector ns by nh ]
%           The derivative of dGds wrt h at the scalar time point t.
%       .d2Gdxdq [ handle @(t,sol) returns real vector nq by nx ]
%           The derivative of dGdq wrt x at the scalar time point t.
%       .d2Gdkdq [ handle @(t,sol) returns real vector nq by nk ]
%           The derivative of dGdq wrt k at the scalar time point t.
%       .d2Gdsdq [ handle @(t,sol) returns real vector nq by ns ]
%           The derivative of dGdq wrt s at the scalar time point t.
%       .d2Gdhdq [ handle @(t,sol) returns real vector nq by nh ]
%           The derivative of dGdq wrt h at the scalar time point t.
%       .d2Gdxdh [ handle @(t,sol) returns real vector nh by nx ]
%           The derivative of dGdh wrt x at the scalar time point t.
%       .d2Gdkdh [ handle @(t,sol) returns real vector nh by nk ]
%           The derivative of dGdh wrt k at the scalar time point t.
%       .d2Gdsdh [ handle @(t,sol) returns real vector nh by ns ]
%           The derivative of dGdh wrt s at the scalar time point t.
%       .d2Gdqdh [ handle @(t,sol) returns real vector nh by nq ]
%           The derivative of dGdh wrt q at the scalar time point t.
%
%       The following are defined for information-theory-based structures.
%       .p [ handle @(sol) returns real scalar ]
%           The likelihood function for this information objective.
%       .logp [ handle @(sol) returns real scalar ]
%           The log likelihood for this information objective.
%       .F [ handle @(sol) 
%            returns symmetric positive semidefinite matrix nT by nT ]
%           Returns the Fisher information matrix for the simulation
%           provided
%       .Fn [ handle @(sol,T) 
%             returns symmetric positive semidefinite matrix nT by nT ]
%           The Fisher information matrix normalized into log parameter
%           space according to active parameters T.
%
%       The following are defined for data objective functions only.
%       .AddData [ handle @(sol) returns objective struct scalar ]
%           Use the simulation to sample data and add it to the objective

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if nargin < 1
    dims = 1;
end

dims = row(dims);

if isscalar(dims)
    dims = [dims, 1];
end

% Objective structure
obj.Type = 'Objective.Data';
obj.Name = 'UnamedObjective';

% Objective function control parameters
obj.Continuous = false; % Continuous g here is 0
obj.Complex = false; % This objective function does not need it
obj.tF = 0; % The latest time at which this objective needs to be run
obj.DiscreteTimes = zeros(1,0); % Times at which discrete measurements are taken

% Event handling
obj.ne = 0;
obj.Events = @empty_event;
obj.IsFinished = @(sol)true;

% Continuous objective function
obj.g = @g;
obj.dgdx = @dgdx;
obj.d2gdx2 = @d2gdx2;

% Discrete objective function
obj.G = @G;
obj.dGdx = @dGdx;
obj.dGdk = @dGdk;
obj.dGds = @dGds;
obj.dGdq = @dGdq;
obj.dGdh = @dGdh;
obj.d2Gdx2 = @d2Gdx2;
obj.d2Gdk2 = @d2Gdk2;
obj.d2Gds2 = @d2Gds2;
obj.d2Gdq2 = @d2Gdq2;
obj.d2Gdh2 = @d2Gdh2;
obj.d2Gdkdx = @d2Gdkdx;
obj.d2Gdsdx = @d2Gdsdx;
obj.d2Gdqdx = @d2Gdqdx;
obj.d2Gdhdx = @d2Gdhdx;
obj.d2Gdxdk = @d2Gdxdk;
obj.d2Gdsdk = @d2Gdsdk;
obj.d2Gdqdk = @d2Gdqdk;
obj.d2Gdhdk = @d2Gdhdk;
obj.d2Gdxds = @d2Gdxds;
obj.d2Gdkds = @d2Gdkds;
obj.d2Gdqds = @d2Gdqds;
obj.d2Gdhds = @d2Gdhds;
obj.d2Gdxdq = @d2Gdxdq;
obj.d2Gdkdq = @d2Gdkdq;
obj.d2Gdsdq = @d2Gdsdq;
obj.d2Gdhdq = @d2Gdhdq;

% Information theory
obj.p = @p;
obj.logp = @logp;
obj.F = @F;
obj.Fn = @Fn;

% Data manipulation
obj.AddData = @AddData;

% Copy the objective structure
obj = repmat(obj, dims);

end

function [value, is_terminal, direction] = empty_event(t, y)
value = zeros(0,1);
is_terminal = false(0,1);
direction = zeros(0,1);
end

function val = g(t,x,u)
val = 0;
end
function val = dgdx(t,x,u)
nx = numel(x);
val = sparse(nx, 1);
end
function val = d2gdx2(t,x,u)
nx = numel(x);
val = sparse(nx, nx);
end

function [val, discreteTimes] = G(sol)
val = 0;
discreteTimes = zeros(0,1);
end
function val = dGdx(t,sol)
nx = size(sol.C1,2);
val = sparse(nx, 1);
end
function val = dGdk(t,sol)
nk = numel(sol.k);
val = sparse(nk, 1);
end
function val = dGds(t,sol)
ns = numel(sol.s);
val = sparse(ns, 1);
end
function val = dGdq(t,sol)
nq = numel(sol.q);
val = sparse(nq, 1);
end
function val = dGdh(t,sol)
nh = numel(sol.h);
val = sparse(nh, 1);
end
function val = d2Gdx2(t,sol)
nx = size(sol.C1,2);
val = sparse(nx, nx);
end
function val = d2Gdk2(t,sol)
nk = numel(sol.k);
val = sparse(nk, nk);
end
function val = d2Gds2(t,sol)
ns = numel(sol.s);
val = sparse(ns, ns);
end
function val = d2Gdq2(t,sol)
nq = numel(sol.q);
val = sparse(nq, nq);
end
function val = d2Gdh2(t,sol)
nh = numel(sol.h);
val = sparse(nh, nh);
end
function val = d2Gdkdx(t,sol)
nk = numel(sol.k);
nx = size(sol.C1,2);
val = sparse(nx, nk);
end
function val = d2Gdsdx(t,sol)
ns = numel(sol.s);
nx = size(sol.C1,2);
val = sparse(nx, ns);
end
function val = d2Gdqdx(t,sol)
nq = numel(sol.q);
nx = size(sol.C1,2);
val = sparse(nx, nq);
end
function val = d2Gdhdx(t,sol)
nh = numel(sol.h);
nx = size(sol.C1,2);
val = sparse(nx, nh);
end
function val = d2Gdxdk(t,sol)
nk = numel(sol.k);
nx = size(sol.C1,2);
val = sparse(nk, nx);
end
function val = d2Gdsdk(t,sol)
ns = numel(sol.s);
nk = numel(sol.k);
val = sparse(nk, ns);
end
function val = d2Gdqdk(t,sol)
nk = numel(sol.k);
nq = numel(sol.q);
val = sparse(nk, nq);
end
function val = d2Gdhdk(t,sol)
nk = numel(sol.k);
nh = numel(sol.h);
val = sparse(nk, nh);
end
function val = d2Gdxds(t,sol)
ns = numel(sol.s);
nx = size(sol.C1,2);
val = sparse(ns, nx);
end
function val = d2Gdkds(t,sol)
nk = numel(sol.k);
ns = numel(sol.s);
val = sparse(ns, nk);
end
function val = d2Gdqds(t,sol)
ns = numel(sol.s);
nq = numel(sol.q);
val = sparse(ns, nq);
end
function val = d2Gdhds(t,sol)
ns = numel(sol.s);
nh = numel(sol.h);
val = sparse(ns, nh);
end
function val = d2Gdxdq(t,sol)
nq = numel(sol.q);
nx = size(sol.C1,2);
val = sparse(nq, nx);
end
function val = d2Gdkdq(t,sol)
nk = numel(sol.k);
nq = numel(sol.q);
val = sparse(nq, nk);
end
function val = d2Gdsdq(t,sol)
ns = numel(sol.s);
nq = numel(sol.q);
val = sparse(nq, ns);
end
function val = d2Gdhdq(t,sol)
ns = numel(sol.s);
nh = numel(sol.h);
val = sparse(nh, ns);
end
function val = d2Gdxdh(t,sol)
nh = numel(sol.h);
nx = size(sol.C1,2);
val = sparse(nh, nx);
end
function val = d2Gdkdh(t,sol)
nk = numel(sol.k);
nh = numel(sol.h);
val = sparse(nh, nk);
end
function val = d2Gdsdh(t,sol)
ns = numel(sol.s);
nh = numel(sol.h);
val = sparse(nh, ns);
end
function val = d2Gdqdh(t,sol)
nq = numel(sol.q);
nh = numel(sol.h);
val = sparse(nh, nq);
end

% For Objective.Information
function val = p(sol)
val = 1;
end
function val = logp(sol)
val = 0;
end
function val = F(sol)
nTk = sum(sol.UseParams);
nTs = sum(sol.UseSeeds);
nTq = sum(sol.UseInputControls);
nTh = sum(sol.UseDoseControls);
nT  = nTk + nTs + nTq + nTh;
val = zeros(nT,nT);
end
function val = Fn(sol)
nTk = sum(sol.UseParams);
nTs = sum(sol.UseSeeds);
nTq = sum(sol.UseInputControls);
nTh = sum(sol.UseDoseControls);
nT  = nTk + nTs + nTq + nTh;
val = zeros(nT,nT);
end

% For Objective.Data
function obj = AddData(sol)
obj = Gzero();
end
