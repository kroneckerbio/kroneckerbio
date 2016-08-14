function obj = objectiveZero(dims)
%objectiveZero Default structure for KroneckerBio objective functions
% 
%   obj = objectiveZero(dims)
%
%   This function returns a objective function that always has a value of
%   0. It is used to fill in the fields with appropriately zero defaults
%   when they are not provided. It also serves as the documentation for the
%   objective function structure
%
%   Inputs
%   dims: [ nonnegative integer vector ]
%       The size of the blank objective structure array
%
%   Outputs
%   obj: [ Kronecker objective structure scalar ]
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
%       .g [ handle @(t,y) returns real scalar ]
%           TODO: Not implemented
%           This is the continuous objective function. This function
%           returns the derivative of the goal function with respect to
%           time. The goal value comes from the numerical integration of
%           this function.
%       .dgdy [ handle @(t,y) returns real vector ny ]
%           The partial derivative of g wrt y
%       .d2gdy2 [ handle @(t,y) returns real matrix ny by ny ]
%           The partial derivative of dgdy wrt y
%       .G [ handle @(int) returns real scalar ]
%           This is the discrete objective function. This function
%           evaluates the integration of the system and returns the
%           objective value. If Complex is true, int will be a dense
%           integration with function handles. Otherwise, int will only
%           provide the solution at DiscreteTimes.
%       .dGdy [ handle @(t,int) returns real vector ny ]
%           The derivative of G wrt y at the scalar time point t.
%       .dGdk [ handle @(int) returns real vector nk ]
%           The derivative of G wrt k at the scalar time point t.
%       .dGds [ handle @(int) returns real vector ns ]
%           The derivative of G wrt s at the scalar time point t.
%       .dGdq [ handle @(t,int) returns real vector nq ]
%           The derivative of G wrt q at the scalar time point t.
%       .dGdh [ handle @(t,int) returns real vector nh ]
%           The derivative of G wrt h at the scalar time point t.
%       .d2Gdy2 [ handle @(t,int) returns real vector ny by ny ]
%           The derivative of dGdy wrt y at the scalar time point t.
%       .d2Gdk2 [ handle @(t,int) returns real vector nk by nk ]
%           The derivative of dGdk wrt k at the scalar time point t.
%       .d2Gds2 [ handle @(t,int) returns real vector ns by ns ]
%           The derivative of dGds wrt s at the scalar time point t.
%       .d2Gdq2 [ handle @(t,int) returns real vector nq by nq ]
%           The derivative of dGdq wrt q at the scalar time point t.
%       .d2Gdh2 [ handle @(t,int) returns real vector nh by nh ]
%           The derivative of dGdq wrt h at the scalar time point t.
%       .d2Gdkdy [ handle @(t,int) returns real vector ny by nk ]
%           The derivative of dGdy wrt k at the scalar time point t.
%       .d2Gdsdy [ handle @(t,int) returns real vector ny by ns ]
%           The derivative of dGdy wrt s at the scalar time point t.
%       .d2Gdqdy [ handle @(t,int) returns real vector ny by nq ]
%           The derivative of dGdy wrt q at the scalar time point t.
%       .d2Gdhdy [ handle @(t,int) returns real vector ny by nh ]
%           The derivative of dGdy wrt h at the scalar time point t.
%       .d2Gdydk [ handle @(t,int) returns real vector nk by ny ]
%           The derivative of dGdk wrt y at the scalar time point t.
%       .d2Gdsdk [ handle @(t,int) returns real vector nk by s ]
%           The derivative of dGdk wrt s at the scalar time point t.
%       .d2Gdqdk [ handle @(t,int) returns real vector nk by nq ]
%           The derivative of dGdk wrt q at the scalar time point t.
%       .d2Gdhdk [ handle @(t,int) returns real vector nk by nh ]
%           The derivative of dGdk wrt h at the scalar time point t.
%       .d2Gdyds [ handle @(t,int) returns real vector ns by ny ]
%           The derivative of dGds wrt y at the scalar time point t.
%       .d2Gdkds [ handle @(t,int) returns real vector ns by nk ]
%           The derivative of dGds wrt k at the scalar time point t.
%       .d2Gdqds [ handle @(t,int) returns real vector ns by ns ]
%           The derivative of dGds wrt q at the scalar time point t.
%       .d2Gdhds [ handle @(t,int) returns real vector ns by nh ]
%           The derivative of dGds wrt h at the scalar time point t.
%       .d2Gdydq [ handle @(t,int) returns real vector nq by ny ]
%           The derivative of dGdq wrt y at the scalar time point t.
%       .d2Gdkdq [ handle @(t,int) returns real vector nq by nk ]
%           The derivative of dGdq wrt k at the scalar time point t.
%       .d2Gdsdq [ handle @(t,int) returns real vector nq by ns ]
%           The derivative of dGdq wrt s at the scalar time point t.
%       .d2Gdhdq [ handle @(t,int) returns real vector nq by nh ]
%           The derivative of dGdq wrt h at the scalar time point t.
%       .d2Gdydh [ handle @(t,int) returns real vector nh by ny ]
%           The derivative of dGdh wrt y at the scalar time point t.
%       .d2Gdkdh [ handle @(t,int) returns real vector nh by nk ]
%           The derivative of dGdh wrt k at the scalar time point t.
%       .d2Gdsdh [ handle @(t,int) returns real vector nh by ns ]
%           The derivative of dGdh wrt s at the scalar time point t.
%       .d2Gdqdh [ handle @(t,int) returns real vector nh by nq ]
%           The derivative of dGdh wrt q at the scalar time point t.
%
%       The following are defined for information-theory-based structures.
%       .p [ handle @(int) returns real scalar ]
%           The likelihood function for this information objective.
%       .logp [ handle @(int) returns real scalar ]
%           The log likelihood for this information objective.
%       .F [ handle @(int) 
%            returns symmetric positive semidefinite matrix nT by nT ]
%           Returns the Fisher information matrix for the integration
%           provided

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if nargin < 1
    dims = 1;
end

dims = row(dims);

if isscalar(dims)
    dims = [dims, 1];
end

% Inherit observation
obj = observationZero();

% Objective structure
obj.Type = 'Objective.Data.Zero';
obj.Name = 'UnamedObjective';

% Objective function control parameters
obj.Continuous = false; % Continuous g here is 0

% Continuous objective function
obj.g = @g;
obj.dgdy = @dgdy;
obj.d2gdy2 = @d2gdy2;

% Discrete objective function
obj.G = @G;
obj.dGdy = @dGdy;
obj.d2Gdy2 = @d2Gdy2;

obj.dGdk = @dGdk;
obj.dGds = @dGds;
obj.dGdq = @dGdq;
obj.dGdh = @dGdh;
obj.d2Gdk2 = @d2Gdk2;
obj.d2Gds2 = @d2Gds2;
obj.d2Gdq2 = @d2Gdq2;
obj.d2Gdh2 = @d2Gdh2;
obj.d2Gdsdk = @d2Gdsdk;
obj.d2Gdqdk = @d2Gdqdk;
obj.d2Gdhdk = @d2Gdhdk;
obj.d2Gdkds = @d2Gdkds;
obj.d2Gdqds = @d2Gdqds;
obj.d2Gdhds = @d2Gdhds;
obj.d2Gdkdq = @d2Gdkdq;
obj.d2Gdsdq = @d2Gdsdq;
obj.d2Gdhdq = @d2Gdhdq;

% lsqnonlin functions
obj.err = @err;
obj.derrdT = @derrdT;

% Information theory
obj.p = @p;
obj.logp = @logp;
obj.F = @F;

% Copy the objective structure
obj = repmat(obj, dims);
end

function val = g(t,y)
val = 0;
end
function val = dgdy(t,y)
ny = numel(y);
val = sparse(ny,1);
end
function val = d2gdy2(t,y)
ny = numel(y);
val = sparse(ny,ny);
end

function [val,discreteTimes] = G(int)
val = 0;
discreteTimes = zeros(0,1);
end
function val = dGdy(t,int)
ny = int.ny;
val = sparse(ny,1);
end
function val = d2Gdy2(t,int)
ny = int.ny;
val = sparse(ny,ny);
end

function val = dGdk(int)
nk = int.nk;
val = sparse(nk,1);
end
function val = dGds(int)
ns = int.ns;
val = sparse(ns,1);
end
function val = dGdq(int)
nq = numel(int.q);
val = sparse(nq,1);
end
function val = dGdh(int)
nh = int.nh;
val = sparse(nh,1);
end
function val = d2Gdk2(int)
nk = int.nk;
val = sparse(nk,nk);
end
function val = d2Gds2(int)
ns = int.ns;
val = sparse(ns,ns);
end
function val = d2Gdq2(int)
nq = numel(int.q);
val = sparse(nq,nq);
end
function val = d2Gdh2(int)
nh = int.nh;
val = sparse(nh,nh);
end
function val = d2Gdsdk(int)
ns = int.ns;
nk = int.nk;
val = sparse(nk,ns);
end
function val = d2Gdqdk(int)
nk = int.nk;
nq = numel(int.q);
val = sparse(nk,nq);
end
function val = d2Gdhdk(int)
nk = int.nk;
nh = int.nh;
val = sparse(nk,nh);
end
function val = d2Gdkds(int)
nk = int.nk;
ns = int.ns;
val = sparse(ns,nk);
end
function val = d2Gdqds(int)
ns = int.ns;
nq = numel(int.q);
val = sparse(ns,nq);
end
function val = d2Gdhds(int)
ns = int.ns;
nh = int.nh;
val = sparse(ns,nh);
end
function val = d2Gdkdq(int)
nk = int.nk;
nq = numel(int.q);
val = sparse(nq,nk);
end
function val = d2Gdsdq(int)
ns = int.ns;
nq = numel(int.q);
val = sparse(nq,ns);
end
function val = d2Gdhdq(int)
ns = int.ns;
nh = int.nh;
val = sparse(nh,ns);
end

function val = err(int)
val = zeros(0,1);
end
function val = derrdT(int)
nT = int.nT;
val = sparse(0,nT);
end

% For Objective.Information
function val = p(int)
val = 1;
end
function val = logp(int)
val = 0;
end
function val = F(int)
nT = int.nT;
val = zeros(nT,nT);
end
