%% Initialize Model
m = InitializeModel('Testing');

%% Add Compartments
% All different dimensions
m = AddCompartment(m, 'v1', 3, '', 1);
m = AddCompartment(m, 'v2', 2, '', 1);

%% Add States
m = AddState(m, 'x1', 'v1');
m = AddState(m, 'x2', 'v1', 10);

%% Add Inputs
m = AddInput(m, 'u1', 'v1', 1);
m = AddInput(m, 'u2', 'v1', @(t,q)2);
m = AddInput(m, 'u3', 'v1', @(t,q)(q), 3);

%% Add Outputs
m = AddOutput(m, 'y1', 'x1', 1);
m = AddOutput(m, 'y2', 'x', 2);
m = AddOutput(m, 'y3', 'u', 3);
m = AddOutput(m, 'y4', '\...$', 4);
m = AddOutput(m, 'y5', '', 5);

%% Add Parameters
m = AddParameter(m, 'k1', 1);
m = AddParameter(m, 'k2', 2);
m = AddParameter(m, 'k3', 3);

%% Add Reactions
m = AddReaction(m, 'r7', '', 'x1', '', '', '', 'k1');
m = AddReaction(m, 'r8', '', 'x1', '', 'x2', '', 'k1');
m = AddReaction(m, 'r9', '', 'x1', '', 'x1', 'x2', 'k1');
m = AddReaction(m, 'r10', '', 'x1', '', 'u3', 'x2', 'k1');
m = AddReaction(m, 'r11', '', 'x1', '', 'u3', 'u1', 'k1');
m = AddReaction(m, 'r12', '', 'x1', '', 'u3', '', 'k1');

%% Finalize Model
m = FinalizeModel(m);

%% Overwrite components
m = AddCompartment(m, 'v2', 2, '', 4);
m = AddState(m, 'x1', 'v1', 14);
m = AddInput(m, 'u2', 'v1', @(t,q)2);
m = AddOutput(m, 'y2', 'x', 1);
m = AddParameter(m, 'k1', 2);
m = AddReaction(m, 'r9', '', 'x1', '', 'x1', 'x2', 'k1');

%% Finalize Model
m = FinalizeModel(m);
