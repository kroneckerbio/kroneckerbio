function [H, All] = FiniteObjectiveHessian(m, con, obj, opts)
%FiniteObjectiveHessian Compute the hessian of an objective function by
%   the finite difference approximation
%
%   Mathematically: H = dG/dp2 or H = dG/dlnp2
%
%   [H, All] = FiniteObjectiveHessian(m, con, obj, opts)
%
%   Inputs:
%       m    - The KroneckerBio model that will be simulated
%       con  - A structure vector of the experimental conditions under
%              which the hessian will be evaluated
%       obj  - A structure array of the objective functions under which the
%              hessian will be evaluated. Each row of obj is matched to the
%              corresponding entry in con
%       opts - Optional function options
%           .useParams   - Vector of indexes indicating the rate constants
%                          whose sensitivities will be considered
%           .useICs      - Vector of indexes indicating the initial
%                          concentrations whose sensitivities will be
%                          considered
%           .UseModelICs - Boolean that determines whether to use the
%                          initial concentrations of the model or the
%                          conditions. Default = true
%           .Normalized  - Boolean that determines if the simple hessian or
%                          the normalized hessian will be computed. The
%                          normalized hessian is normalized with respect to
%                          the values of the parameters. Default = true
%           .ImaginaryStep - Scalar boolean. If set to true, use an
%                           imaginary finite difference step. This has the
%                           advantage of not requiring the subtraction of
%                           two large numbers, increasing the stability of
%                           the result, but comes at the cost of larger
%                           computational cost from evaluating expressions
%                           with complex numbers. If set to false (the
%                           default), a real finite difference step is
%                           used.
%           .Verbose     - Print progress to command window
%   Outputs:
%       H = FiniteObjectiveHessian(m, con, obj, ...)
%           H - The hessian as an array
%
%       [H, All] = FiniteObjectiveHessian(m, con, obj, ...)
%           All - The hessian is the sum of hessians from all experiements,
%           but All returns the hessians in a size(obj) cell array
%
%   Additional info:
%   - The experimental condition vector can also be a cell vector
%   - The objective function array can also be a cell array. Empty entries 
%   in the cell array and entries in the structure array with empty
%   values are ignored. This way, conditions can have different numbers of
%   objective functions associated with them.

% (c) 2015 David R Hagen and Bruce Tidor
% This software is released under the GNU GPLv2 or later.

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
end

assert(nargin >= 3, 'KroneckerBio:ObjectiveGradient:TooFewInputs', 'ObjectiveGradient requires at least 3 input arguments')
assert(isscalar(m), 'KroneckerBio:ObjectiveGradient:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = [];
defaultOpts.AbsTol         = [];

defaultOpts.ImaginaryStep  = false;

defaultOpts.UseParams        = 1:m.nk;
defaultOpts.UseSeeds         = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls  = [];

defaultOpts.ObjWeights     = ones(size(obj));

defaultOpts.Normalized     = true;
defaultOpts.UseAdjoint     = false;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
ns = m.ns;
nk = m.nk;
n_con = numel(con);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, n_con);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, n_con, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, n_con, cat(1,con.nh));

nT = nTk + nTs + nTq + nTh;

% Store starting parameter sets
T0 = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

% Refresh conditions and objectives
con = refreshCon(m, con);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, n_con, opts.UseAdjoint, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

%% Loop through conditions
H = zeros(nT,nT);

% Initial value
if verbose; fprintf('Initial round\n'); end
[unused, D] = computeObjGrad(m, con, obj, opts);

for iT = 1:nT
    if verbose; fprintf('Step %d of %d\n', iT, nT); end
    
    % Set baseline parameters
    T_i = T0(iT);
    T_up = T0;
    
    % Change current parameter by finite amount
    if opts.ImaginaryStep
        imagfactor = 1i;
    else
        imagfactor = 1;
    end
    if opts.Normalized
        normfactor = T_i;
    else
        normfactor = 1;
    end
    stepsize = 1e-8;
    diff = normfactor * imagfactor * stepsize;
    
    % Compute objective values
    T_up(iT) = T_up(iT) + diff;
    [m, con] = updateAll(m, con, T_up, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
    [unused, D_up] = computeObjGrad(m, con, obj, opts);

    % Compute D
    if opts.ImaginaryStep
        H(:,iT) = imag(D_up) ./ imag(diff);
    else
        H(:,iT) = (D_up - D) ./ diff;
    end
    if opts.Normalized
        H(:,iT) = T_i * T0 .* H(:,iT);
    end
    
end
