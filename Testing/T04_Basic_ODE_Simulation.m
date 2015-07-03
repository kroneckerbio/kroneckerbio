%% Construct equilibrium experiment A + B <-> C
m = InitializeModelMassActionAmount('Equilibrium');

m = AddCompartment(m, 'Solution', 3, 1);

m = AddState(m, 'A', 'Solution', 1);
m = AddState(m, 'B', 'Solution', 2);
m = AddState(m, 'C', 'Solution', 0);

m = AddOutput(m, 'A', 'A');
m = AddOutput(m, 'B', 'B');
m = AddOutput(m, 'C', 'C');

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);

m = AddReaction(m, '', 'A', 'B', 'C', '', 'kf', 'kr');

m = FinalizeModel(m);

%% Construct experiment
con = experimentInitialValue(m, [], [], [], 'InitialValueExperiment');

%% Simulate
tF = 1; % final time
sim1 = SimulateSystem(m, con, tF);

figure
times = linspace(0, tF, 100);
plot(times, sim1.x(times))
legend('A','B','C')
xlabel('Time')
ylabel('Amount')

%% Simulate Sensitivities
sim2 = SimulateSensitivity(m, con, tF);

%% Simulate Curvature
opts.UseSeeds = [];
sim3 = SimulateCurvature(m, con, tF, opts);

%% Simulate Linear Noise Approximation
% sim4 = SimulateLna(m, con);