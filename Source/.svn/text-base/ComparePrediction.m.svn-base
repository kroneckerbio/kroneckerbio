function result = ComparePrediction(m1, m2, con, com, opts)
%COMPAREPREDICTION Compares the predictions of two models under a set of
%   experimental conditions according to a supplied compare function.
%
%   result = ComparePrediction(m1, m2, con, com, opts)
%
%   Inputs
%       m1   - Kronecker Bio Model
%       m2   - Kronecker Bio Model
%       con  - Vector of Kronecker Bio experimental conditions
%       com  - Scalar Kronecker Bio compare structure
%       opts - Scalar structure of options
%           UseModelICs - Logical scalar indicating that the model ICs
%                         should be use in place of the con ICs.
%                         Default = true
%           AbsTol - Scalar absolute integration tolerance for the system. 
%                    Default = 1e-6
%           RelTol - Scalar relative integration tolerance for the system.
%                    Default = 1e-4
%           PredictionType - Passed to com, which can modulate the results
%                            it returns accordingly.
%                            Default = []
%           Verbose - Scalar depth of verbosity. Higher numbers lead to
%                     more detailed print out.
%
%   Outputs
%       results - A matrix length of con by length of whatever com returns
%
%   What this function does is largely up to the com structure. It
%   simulates both models under each experimental condition and queries com
%   to process the simulation results. Each row in result contains the
%   results from a single experimental condition.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 5
    opts = [];
end

assert(isscalar(m1), 'KroneckerBio:FitObjective:MoreThanOneModel', 'The model structure must be scalar')
assert(isscalar(m2), 'KroneckerBio:FitObjective:MoreThanOneModel', 'The model structure must be scalar')
assert(isscalar(com), 'KroneckerBio:FitObjective:MoreThanOneModel', 'The comparison structure must be scalar')

%% Options
% Default options
defaultOpts.AbsTol          = NaN;
defaultOpts.RelTol          = NaN;
defaultOpts.UseModelICs     = true;
defaultOpts.PredictionType  = [];
defaultOpts.Verbose         = 1;

opts = mergestruct(defaultOpts, opts);
verbose = logical(opts.Verbose);

% Constants
nX1 = m1.nX;
nX2 = m2.nX;
nCon = numel(con);

%% Integration type: simple, continuous, complex, or both
% Fix integration type
opts.tGet = com.discreteTimes;
opts.complex = com.complex;

%% Tolerances
if isempty(opts.RelTol) || isnan(opts.RelTol)
    opts.RelTol = 1e-6;
end

% Fix AbsTol to be a cell array of vectors appropriate to the problem
if ~iscell(opts.AbsTol)
    abstol1 = fixAbsTol(opts.AbsTol, 1, zeros(nCon,1), nX1, nCon);
    abstol2 = fixAbsTol(opts.AbsTol, 1, zeros(nCon,1), nX2, nCon);
end

%% Options structure for integration
opts1 = opts;
opts2 = opts;
opts1.AbsTol = abstol1;
opts2.AbsTol = abstol2;
intOpts1 = opts1;
intOpts2 = opts2;

%% Loop over all experiments
result = com.empty(nCon);
for iCon = 1:nCon;
    if verbose; fprintf(['Comparing predictions for ' con(iCon).name '...']); end
    
    % Modify opts structure
    intOpts1.AbsTol = opts1.AbsTol{iCon};
    intOpts2.AbsTol = opts2.AbsTol{iCon};
    
    % Simulate both systems
    if opts.complex
        sol1 = integrateSys(m1, con(iCon), intOpts1);
        sol2 = integrateSys(m2, con(iCon), intOpts2);
    else
        sol1 = integrateSysDisc(m1, con(iCon), intOpts1);
        sol2 = integrateSysDisc(m2, con(iCon), intOpts2);
    end
    
    % Append the appropriate output vector
    sol1.c = m1.c;
    sol2.c = m2.c;
    
    result(iCon,:) = com.compare(sol1, sol2, opts.PredictionType);
end
