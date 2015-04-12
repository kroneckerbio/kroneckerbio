function absTol = GoodAbsTol(m, con, sd, opts)
%GoodAbsTol Make a reasonable estimate as what the absolute tolerance on
%   the species and sensitivities for objective calculations
%
%   absTolRatio = GoodAbsTol(m, con, sd, opts)
%
%   This function uses the lowest uncertainty in the outputs, as given by
%   sd at an output value of 0, as the floor on the relavent output value.
%   This output value is then propogated to a floor on the relavent species
%   values. This floor is the AbsTolRatio. AbsTolRatio times RelTol is a
%   good guess as to what should be used as AbsTol for ODE integration of
%   the states and sensitivities.
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model
%   con: [ experiment struct vector ]
%       The experimental conditions 
%   sd: [ handle @(t,yInd,yVal) returns positive scalar ]
%       The standard deviation as a function of the time, output index, and
%       output value.
%   opts: [ options struct scalar ]
%       Optional
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Which kinetic parameters are active
%       .UseSeeds [ logical matrix nx by nCon | logical vector nx |
%                   positive integer vector {[]} ]
%           Which seed parameters are active
%       .UseInputControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Which input control parameters arae active
%       .UseDoseControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Which dose control parameters are active
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration. This is used to scale
%           the good absolute tolerances.
%
%   Outputs
%   absTol: [ struct vector nCon ]
%       An AbsTol structure that is valid for all simulations

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Constants
nk = m.nk;
ns = m.ns;
nx = m.nx;
ny = m.ny;
nCon = numel(con);

% Default options
defaultOpts.RelTol         = [];

defaultOpts.UseParams        = 1:m.nk;
defaultOpts.UseSeeds         = nan;
defaultOpts.UseInputControls = nan;
defaultOpts.UseDoseControls  = nan;

opts = mergestruct(defaultOpts, opts);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseICs is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, nCon);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, nCon, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, nCon, cat(1,con.nh));

nT = nTk + nTs + nTq + nTh;

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

%% xGoodRatio
% Get the typical uncertainties for the outputs
yFloor = zeros(ny,1);    % The lowest uncertainty for this output
yStandard = zeros(ny,1); % A typical uncertainty for this output
for iy = 1:ny
    yFloor(iy) = sd(0, iy, 0);
    yStandard(iy) = sd(0, iy, yFloor(iy)*10) ./ yFloor(iy)*10;
end

% A good abstol ratio is the ratio of the floor to the standard uncertainty
yRatio = yFloor ./ yStandard;

% Transfer the floor back to the species
% Out of all the outputs that a species affects, we should take the most
% conservative floor.
xRatio = bsxfun(@rdivide, yRatio, m.dydx(0, ones(m.nx,1), ones(m.nu,1))); % Propogation of floors
xRatio = min(xRatio, [], 1); % Keep most conservative
xRatio = vec(xRatio); % Straighten it

% Assume that the effect between species is one-to-one and can travel an
% infinite distance along the network.
% net = NetworkDistance(m); % Distance between species
% net = spones(net); % Convert to one-to-one mapping
% xGoodRatio = bsxfun(@rdivide, xFloor, net); % Propogation
xGoodRatio = repmat(min(xRatio, [], 1), nx,1); % Keep most conservative and assign it to all
xGoodRatio = vec(xGoodRatio); % Straighten it

% The objective is about 1 for a single data point
objectiveGoodRatio = 1;

%% Assemble the vectors for each experiment
absTol = emptystruct(nCon, 'System', 'Objective', 'Sensitivity', 'Adjoint', 'Gradient', 'Doublesensitivity');
for iCon = 1:nCon
    % Experiment specific parameters
    T = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
    
    % dxdTGoodRatio
    % If the model is linear and x >= 0 for all T >= 0, then the
    % maximum value of dx/dT is x/T. A good ratio would therefore be
    % x/T
    dxdTGoodRatio = bsxfun(@rdivide, xGoodRatio, T.'); % x_T
    dxdTGoodRatio = vec(dxdTGoodRatio); % xT_
    
    % d2xdT2GoodRatio
    % If the model is quadratic and x >= 0 for all T >= 0, then the
    % maximum value of d2x/dT2 is 2*x/T^2. A good ratio would therefore
    % be 2*x/T^2
    d2xdT2GoodRatio = bsxfun(@rdivide, xGoodRatio, bsxfun(@times, T.', reshape(T, 1,1,nT))); % x_T_T
    d2xdT2GoodRatio = vec(d2xdT2GoodRatio); % xTT_
    
    % D ratio
    % These elements are not computed in log space so their absolute
    % value will be inversely proportational to theta.
    DGoodRatio = 1 ./ T;
    
    % lambda ratio
    % I have no idea what to set this to. Using the other variables to
    % control the error seems to be effective
    lAbsTol = inf(nx,1);
    
    % Distribute heuristic AbsTols
    absTol(iCon).System = opts.RelTol * xGoodRatio;
    absTol(iCon).Objective = opts.RelTol * [xGoodRatio; objectiveGoodRatio];
    absTol(iCon).Sensitivity = opts.RelTol * [xGoodRatio; dxdTGoodRatio];
    absTol(iCon).Gradient = opts.RelTol * [xGoodRatio; objectiveGoodRatio; dxdTGoodRatio; DGoodRatio];
    absTol(iCon).Adjoint = opts.RelTol * [xGoodRatio; lAbsTol; DGoodRatio];
    absTol(iCon).Curvature = opts.RelTol * [xGoodRatio; dxdTGoodRatio; d2xdT2GoodRatio];
end
