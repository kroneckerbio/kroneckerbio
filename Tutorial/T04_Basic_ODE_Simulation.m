%% T04 Basic ODE Simulation
% Run a simple simulation and plot the results.
% Do some additional analysis: sensitivity and curvature calculation, and linear
%   noise approximation to make sure they work.

%% Load equilibrium experiment A + B <-> C
current_path = fileparts(mfilename('fullpath'));
m = LoadModelMassAction(fullfile(current_path, '../Testing/Equilibrium.txt'));

%% Construct experiment
con = experimentInitialValue(m, [], [], [], 'InitialValueExperiment');

%% Simulate
tF = 1; % final time
sim1 = SimulateSystem(m, con, tF);

% Plot result
figure
times = linspace(0, tF, 100);
plot(times, sim1.x(times))
legend('A','B','C')
xlabel('Time')
ylabel('Amount')

%% Simulate Sensitivities
sim2 = SimulateSensitivity(m, con, tF);

%% Simulate Curvature
sim3 = SimulateCurvature(m, con, tF);

%% Simulate Linear Noise Approximation
obs = observationAll(tF);
sim4 = SimulateLna(m, con, obs);
