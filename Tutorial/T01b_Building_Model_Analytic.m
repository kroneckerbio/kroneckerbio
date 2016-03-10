%% T01b Building an Analytic Model
% Build an analytic model, demonstrating different ways of calling the
%   component-adding functions.
% Analytic models are more general than massaction models, allowing rates as
%   arbitrary functions of states, inputs, and parameters; outputs as arbitray
%   functions of states, inputs, and parameters; and initial conditions as
%   arbitrary functions of seeds.

%% Initialize Model
m = InitializeModelAnalytic('TestAnalyticModel');

%% Add Compartments
m = AddCompartment(m, 'v1', 3, 1);
m = AddCompartment(m, 'v2', 3, 1);
m = AddCompartment(m, 'v3', 2, 1);
m = AddCompartment(m, 'v4', 1, 1);
m = AddCompartment(m, 'v5', 0, 1);

%% Add Seeds
m = AddSeed(m, 's1', 5);
m = AddSeed(m, 's2', 3);
m = AddSeed(m, 's3:s4', 4);

%% Add States
% Demonstrate initial conditions as arbitrary expressions
m = AddState(m, 'x0', 'v1', 1);
m = AddState(m, 'x0', 'v2', '2*s1 + 1.5^("s3:s4")');
m = AddState(m, 'x1', 'v1');
m = AddState(m, 'x2', 'v1', 10);
m = AddState(m, 'x3', 'v1', 0);
m = AddState(m, 'x4', 'v1', 's1');
m = AddState(m, 'x5', 'v1', '2*s1 + 1.5^(s2)');

%% Add Inputs
m = AddInput(m, 'u1', 'v1', 1);
m = AddInput(m, 'u2', 'v1', 2);
m = AddInput(m, 'u3', 'v1', 3);
m = AddInput(m, 'u4', 'v1', 4);

%% Add Parameters
m = AddParameter(m, 'k1', 1);
m = AddParameter(m, 'k2', 2);
m = AddParameter(m, 'k3', 3);
m = AddParameter(m, 'k4', 4);

%% Add Outputs
% Demonstrate outputs as arbitrary expressions
m = AddOutput(m, 'y1', 'x1');
m = AddOutput(m, 'y2', '2*x1 + 3.5*x2');
m = AddOutput(m, 'y3', 'k1*x1');
m = AddOutput(m, 'y4', 'exp(k2*x2)');
m = AddOutput(m, 'y5', 'v1.x0');
m = AddOutput(m, 'y6', 'v2.x0');

%% Add Reactions
% Demonstrate rates as arbitrary expressions
m = AddReaction(m, 'r00', 'x1', {'x3', 'x4'}, 'k1*x5');
m = AddReaction(m, 'r01', 'x1', {'x3', 'x4'}, 'k1 + k2');
m = AddReaction(m, 'r02', 'x1', {'x3', 'x4'}, 'k1/(k2^2 + x1^2)');
m = AddReaction(m, 'r03', 'x1', {'x3', 'x4'}, 'k1*x0/(k2^2 + u2^2)', '', 'v1');
m = AddReaction(m, 'r04', 'x1', {'x3', 'x4'}, 'k1*v1.x0/(k2^2 + u2^2)');
m = AddReaction(m, 'r05', 'x1', {'x3', 'x4'}, 'k1*v1.x0/(k2^2 + u2^2)', '', 'v1');
m = AddReaction(m, 'r06', 'x1', {'x3', 'x4'}, 'exp(k1 + 1.4*u3)');
m = AddReaction(m, 'r1', {}, {}, 'k1');
m = AddReaction(m, 'r2', {}, 'x3', 'k1');
m = AddReaction(m, 'r3', {}, {'x3', 'x4'}, 'k1');
m = AddReaction(m, 'r4', {}, {'u3', 'x4'}, 'k1');
m = AddReaction(m, 'r5', {}, {'u3', 'u4'}, 'k1');
m = AddReaction(m, 'r6', {}, 'u3', 'k1');

m = AddReaction(m, 'r7', 'x1', {}, 'k1');
m = AddReaction(m, 'r8', 'x1', 'x3', 'k1');
m = AddReaction(m, 'r9', 'x1', {'x3', 'x4'}, 'k1');
m = AddReaction(m, 'r10', 'x1', {'u3', 'x4'}, 'k1');
m = AddReaction(m, 'r11', 'x1', {'u3', 'u4'}, 'k1');
m = AddReaction(m, 'r12', 'x1', 'u3', 'k1');

m = AddReaction(m, 'r13', {'x1', 'x2'}, {}, 'k1');
m = AddReaction(m, 'r14', {'x1', 'x2'}, 'x3', 'k1');
m = AddReaction(m, 'r15', {'x1', 'x2'}, {'x3', 'x4'}, 'k1');
m = AddReaction(m, 'r16', {'x1', 'x2'}, {'u3', 'x4'}, 'k1');
m = AddReaction(m, 'r17', {'x1', 'x2'}, {'u3', 'u4'}, 'k1');
m = AddReaction(m, 'r18', {'x1', 'x2'}, 'u3', 'k1');

m = AddReaction(m, 'r19', {'x1', 'u2'}, {}, 'k1');
m = AddReaction(m, 'r20', {'x1', 'u2'}, 'x3', 'k1');
m = AddReaction(m, 'r21', {'x1', 'u2'}, {'x3', 'x4'}, 'k1');
m = AddReaction(m, 'r22', {'x1', 'u2'}, {'u3', 'x4'}, 'k1');
m = AddReaction(m, 'r23', {'x1', 'u2'}, {'u3', 'u4'}, 'k1');
m = AddReaction(m, 'r24', {'x1', 'u2'}, 'u3', 'k1');

m = AddReaction(m, 'r25', {'u1', 'u2'}, {}, 'k1');
m = AddReaction(m, 'r26', {'u1', 'u2'}, 'x3', 'k1');
m = AddReaction(m, 'r27', {'u1', 'u2'}, {'x3', 'x4'}, 'k1');
m = AddReaction(m, 'r28', {'u1', 'u2'}, {'u3', 'x4'}, 'k1');
m = AddReaction(m, 'r29', {'u1', 'u2'}, {'u3', 'u4'}, 'k1');
m = AddReaction(m, 'r30', {'u1', 'u2'}, 'u3', 'k1');

m = AddReaction(m, 'r31', 'u2', {}, 'k1');
m = AddReaction(m, 'r32', 'u2', 'x3', 'k1');
m = AddReaction(m, 'r33', 'u2', {'x3', 'x4'}, 'k1');
m = AddReaction(m, 'r34', 'u2', {'u3', 'x4'}, 'k1');
m = AddReaction(m, 'r35', 'u2', {'u3', 'u4'}, 'k1');
m = AddReaction(m, 'r36', 'u2', 'u3', 'k1');

% Demonstrates reaction compartment
m = AddReaction(m, 'r37', 'x1', {'x1', 'x1'}, 'k1*x1', '', 'v1');

% Demonstrate reaction with > 2 reactants/products
m = AddReaction(m, 'r38', {'x1','u1'}, {'x2','x3','x4'}, 'k1');

% Currently, analytic models don't warn and ignore identical reactions
m = AddReaction(m, 'r22', {'x1', 'u2'}, {'x4', 'u3'}, 'k1');

%% Finalize Model
opts = [];
opts.Verbose = 2;
m = FinalizeModel(m, opts);
