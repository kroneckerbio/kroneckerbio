function m = InitializeModelAnalytic(name)
% Create a symbolic model to store analytic model precursors. Call FinalizeModel
% to covert to final analytic model expected by kroneckerbio.
% Inputs:
%   name [ string {''} ]
% Outputs:
%   m [ symbolic model struct ]

if nargin < 1
    name = [];
end
if isempty(name)
    name = '';
end

m.Type = 'Model.Analytic';
m.Name = fixModelName(name);

m.Compartments = growCompartmentsAnalytic;
m.Parameters   = growParametersAnalytic;
m.Seeds        = growSeedsAnalytic;
m.Inputs       = growInputsAnalytic;
m.States       = growStatesAnalytic;
m.Reactions    = growReactionsAnalytic;
m.Outputs      = growOutputsAnalytic;
m.Rules        = growRulesAnalytic;

m = initializeModelBase(m);

m.add.Compartments = growCompartmentsAnalytic;
m.add.Parameters   = growParametersAnalytic;
m.add.Seeds        = growSeedsAnalytic;
m.add.Inputs       = growInputsAnalytic;
m.add.States       = growStatesAnalytic;
m.add.Reactions    = growReactionsAnalytic;
m.add.Outputs      = growOutputsAnalytic;
m.add.Rules        = growRulesAnalytic;
