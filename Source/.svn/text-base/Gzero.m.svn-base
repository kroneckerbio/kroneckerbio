function obj = Gzero(m)
%Gzero Default structure for Kronecker Bio objective functions
% 
%   obj = Gzero(m)
%
%   This function is a blank objective function. It is used to fill in the
%   fields with appropriately zero defaults when they are not provided.
%
%   Inputs
%   m: [ Kronecker model structure scalar ]
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
%       .Linked [ nonnegative integer scalar ]
%           If this is zero, it means that this objective function acts
%           indepently of any other objective function. That is, it only
%           needs the simulation from one experiment. If this is any other
%           number, then the simulations from all the experiments which
%           have an objective function with this integer attached will be
%           passed as a vector to all those objective functions.
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
%       .dgdk [ handle @(t,x,u) returns real vector nx ]
%           The partial derivative of g wrt k
%       .d2gdx2 [ handle @(t,x,u) returns real matrix nx by nx ]
%           The partial derivative of dgdx wrt x
%       .d2gdk2 [ handle @(t,x,u) returns real matrix nk by nk ]
%           The partial derivative of dgdk wrt k
%       .d2gdxdk [ handle @(t,x,u) returns real matrix nk by nx ]
%           The partial derivative of dgdk wrt x
%       .d2gdkdx [ handle @(t,x,u) returns real matrix nx by nk ]
%           The partial derivative of dgdx wrt k
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
%       .d2Gdx2 [ handle @(t,sol) returns real vector nx by nx ]
%           The derivative of dGdx wrt x at the scalar time point t.
%       .d2Gdk2 [ handle @(t,sol) returns real vector nk by nk ]
%           The derivative of dGdk wrt k at the scalar time point t.
%       .d2Gdxdk [ handle @(t,sol) returns real vector nk by nx ]
%           The derivative of dGdk wrt x at the scalar time point t.
%       .d2Gdkdx [ handle @(t,sol) returns real vector nx by nk ]
%           The derivative of dGdx wrt k at the scalar time point t.
%       .Update [ handle @(m,con,UseParams,UseICs,UseControls)
%                 returns objective struct scalar ]
%           This is a how to modify objective structures that depend on the
%           parameters of the model and experiment.
%
%       The following are defined for information-theory-based structures.
%       .F [ handle @(sol) 
%            returns symmetric positive semidefinite matrix nT by nT ]
%           Returns the Fisher information matrix for the simulation
%           provided
%       .Fn [ handle @(sol,T) 
%             returns symmetric positive semidefinite matrix nT by nT ]
%           The Fisher information matrix normalized into log parameter
%           space according to active parameters T.
%       .p [ handle @(sol) returns real scalar ]
%           The likelihood function for this information objective.
%       .pvalue [ handle @(sol,obj) returns real scalar 0 < pvalue < 1 ]
%           The pvalue for the likelihood function. The highly nonlinear
%           nature of the pvalue function requires that all objective
%           functions of interest be compatiable and evaluated together.
%       .n [ nonegative integer scalar ]
%           The number of degrees of freedom in this objective.
%
%       The following are defined for data objective functions only.
%       .AddData [ handle @(sol) returns objective struct scalar ]
%           Use the simulation to sample noisy data and add it to the
%           objective
%       .AddExpectedData [ handle @(sol) returns objective struct scalar ]
%           Use the simulation to sample data, but without noise, and it to
%           the objective
%       .QuantifyPrediction [ handle @(sol1,sol2,type) returns 
%                             cell vector ]
%           Compares the predicted simulation sol1 with the correct
%           simulation sol2 according to type, which is defined differently
%           for each type of objective.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Special case: return empty structure array if inputs are numeric
if isnumeric(m)
    obj = emptystruct(m, 'Type', 'Name', 'Continuous', 'Complex', 'Linked', 'DiscreteTimes', 'g', 'dgdx', 'dgdk', 'd2gdx2', 'd2gdk2', 'd2gdxdk', 'd2gdkdx', 'G', 'dGdx', 'dGdk', 'd2Gdx2', 'd2Gdk2', 'd2Gdkdx', 'd2Gdxdk', 'F', 'Fn', 'dFndT', 'p', 'pvalue', 'n', 'AddData', 'AddExpectedData', 'QuantifyPrediction', 'Update');
    if ~isempty(obj)
        [obj.Update] = deal(@Update);
    end
    return
end

% Constants
nx = m.nx;
nk = m.nk;

% Gut m
temp = m;
clear m
m.nx = temp.nx;
m.nk = temp.nk;
clear temp

% Objective structure
obj.Type = 'Objective.Data';
obj.Name = 'UnamedObjective';

% Objective function control parameters
obj.Continuous = false; % Continuous g here is 0
obj.Complex = false; % This objective function does not need it
obj.Linked = 0; % This objective function requires only one experiment
obj.DiscreteTimes = zeros(0,1); % Times at which discrete measurements are taken

% Continuous objective function
obj.g        = @g;
obj.dgdx     = @dgdx;
obj.dgdk     = @dgdk;
obj.d2gdx2   = @d2gdx2;
obj.d2gdk2   = @d2gdk2;
obj.d2gdxdk  = @d2gdxdk;
obj.d2gdkdx  = @d2gdkdx;

% Discrete objective function
obj.G       = @G;
obj.dGdx    = @dGdx;
obj.dGdk    = @dGdk;
obj.d2Gdx2  = @d2Gdx2;
obj.d2Gdk2  = @d2Gdk2;
obj.d2Gdkdx = @d2Gdkdx;
obj.d2Gdxdk = @d2Gdxdk;

% Information theory
obj.F  = @F;
obj.Fn = @Fn;
obj.dFndT = @dFndT;
obj.p  = @p;
obj.pvalue = @pvalue;

% Data manipulation
obj.n = 0; % Number of discrete measurements
obj.AddData = @AddData;
obj.AddExpectedData = @AddExpectedData;
obj.QuantifyPrediction = @QuantifyPrediction;

% Update
obj.Update = @Update;

    function val = g(t,x,u)
        val = 0;
    end
    function val = dgdx(t,x,u)
        val = sparse(nx, 1);
    end
    function val = dgdk(t,x,u)
        val = sparse(nk, 1);
    end
    function val = d2gdx2(t,x,u)
        val = sparse(nx, nx);
    end
    function val = d2gdk2(t,x,u)
        val = sparse(nk, nk);
    end
    function val = d2gdkdx(t,x,u)
        val = sparse(nx, nk);
    end
    function val = d2gdxdk(t,x,u)
        val = sparse(nk, nx);
    end

    function [val discreteTimes] = G(sol)
        val = 0;
        discreteTimes = zeros(0,1);
    end
    function val = dGdx(t,sol)
        val = sparse(nx, 1);
    end
    function val = dGdk(t,sol)
        val = sparse(nk, 1);
    end
    function val = d2Gdx2(t,sol)
        val = sparse(nx, nx);
    end
    function val = d2Gdk2(t,sol)
        val = sparse(nk, nk);
    end
    function val = d2Gdkdx(t,sol)
        val = sparse(nx, nk);
    end
    function val = d2Gdxdk(t,sol)
        val = sparse(nk, nx);
    end
    
    % For Objective.Information
    function val = F(dxdTSol)
        nT = (size(dxdTSol.y, 1) - nx) / nx;
        val = zeros(nT,nT);
    end
    function val = Fn(dxdTSol, T)
        nT = (size(dxdTSol.y, 1) - nx) / nx;
        val = zeros(nT,nT);
    end
    function val = dFndT(d2xdT2Sol, T)
        nT = numel(T);
        val = zeros(nT*nT,nT);
    end
    function val = p(sol)
        val = 1;
    end
    function val = pvalue(sol)
        val = 1;
    end

    % For Objective.Data
    function obj = AddData(sol)
        obj = Gzero(m);
    end
    function obj = AddExpectedData(sol)
        obj = Gzero(m);
    end
    function val = QuantifyPrediction(sol1, sol2, type)
        val = [];
    end

    function obj = Update(m,con,useParams,useICs,useControls)
        obj = Gzero(m);
    end

end