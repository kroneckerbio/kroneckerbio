%% Construct equilibrium experiment A + B <-> C
m = InitializeModel('Equilibrium');

m = AddCompartment(m, 'Solution', 3, 1);

m = AddState(m, 'A', 'Solution', 1);
m = AddState(m, 'B', 'Solution', 2);
m = AddState(m, 'C', 'Solution', 0);

m = AddOutput(m, 'A', 'A');
m = AddOutput(m, 'B', 'B');
m = AddOutput(m, 'C', 'C');

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);

m = AddReaction(m, '', '', 'A', 'B', 'C', '', 'kf', 'kr');

m = FinalizeModel(m);

%% Construct experiment
con = Experiment(m, 1, m.s, false, false, [], [], [], [], 'EquilibriumExperiment');

%% Construct objective
outputlist = [1;1;2;2;3;3];
timeslist = [0.1;1;0.1;1;0.1;1];
measurements = [0.6; 0.4; 1.5; 1.3; 0.4; 0.6];
sd = sdLinear(0.05, 0.1);

obj = objectiveWeightedSumOfSquaresNonNeg(outputlist, timeslist, sd, measurements);

%% Fit
mfit = FitObjective(m, con, obj);

%% Linearized parameter uncertainty
F = ObjectiveInformation(mfit, con, obj);
