%% Construct equilibrium experiment A + B <-> C
m = InitializeModel('Equilibrium');

m = AddCompartment(m, 'Solution', 3, '', 1);

m = AddState(m, 'A', 'Solution', 1);
m = AddState(m, 'B', 'Solution', 2);
m = AddState(m, 'C', 'Solution', 0);

m = AddOutput(m, 'A', 'A', 1);
m = AddOutput(m, 'B', 'B', 1);
m = AddOutput(m, 'C', 'C', 1);

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);

m = AddReaction(m, '', '', 'A', 'B', 'C', '', 'kf', 'kr');

m = FinalizeModel(m);

%% Construct experiment
con = Experiment(m, 1, m.x0, false, false, [], [], [], [], 'EquilibriumExperiment');

%% Simulate
sim1 = Simulate(m, con);

%% Simulate Sensitivities
sim2 = SimulateSensitivity(m, con);
