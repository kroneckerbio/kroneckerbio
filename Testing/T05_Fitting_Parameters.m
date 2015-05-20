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
con = experimentInitialValue(m, [], [], [], 'InitialValueExperiment');

%% Construct objective
outputList = [1;1;2;2;3;3];
timesList = [0.1;1;0.1;1;0.1;1];
measurements = [0.6; 0.4; 1.5; 1.3; 0.4; 0.6];
sd = sdLinear(0.05, 0.1);

obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation');
obj = obs.Objective(measurements);

%% Fit
mFit = FitObjective(m, con, obj);

% Display fit results
tF = 1;
times = linspace(0, tF, 100);
simOriginal = SimulateSystem(m, con, tF);
simFit = SimulateSystem(mFit, con, tF);

figure
hold on
plot(times, simOriginal.x(times))
ax = gca;
ax.ColorOrderIndex = 1;
plot(timesList(1:2), reshape(measurements,2,3), '+')
ax = gca;
ax.ColorOrderIndex = 1;
plot(times, simFit.x(times), ':')
hold off
legend('A','B','C','A data','B data','C data','A fit','B fit','C Fit')
xlabel('Time')
ylabel('Amount')

%% Linearized parameter uncertainty
F = ObjectiveInformation(mFit, con, obj);
