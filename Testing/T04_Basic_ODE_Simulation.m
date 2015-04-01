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
con = InitialValueExperiment(m, 1, [], [], [], 'InitialValueExperiment');

%% Simulate
sim1 = Simulate(m, con);

%% Simulate Sensitivities
sim2 = SimulateSensitivity(m, con);

%% Simulate Curvature
sim3 = SimulateCurvature(m, con);

%% Simulate Linear Noise Approximation
sim4 = SimulateLna(m, con);
