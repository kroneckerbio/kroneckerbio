%% Initialize Model
m = InitializeModel('Testing');

%% Add Compartments
% All different dimensions
m = AddCompartment(m, 'v1', 3, '', 1);
m = AddCompartment(m, 'v2', 2, '', 1);
m = AddCompartment(m, 'v3', 1, '', 1);
m = AddCompartment(m, 'v4', 0, '', 0);

% All different expressions
m = AddCompartment(m, 'v5', 3, 'x1', 1);
m = AddCompartment(m, 'v6', 3, 'x', 1);
m = AddCompartment(m, 'v7', 3, 'u', 1);
m = AddCompartment(m, 'v8', 3, '^v8\.', 1);

%% Add States
m = AddState(m, 'x1', 'v1');
m = AddState(m, 'x2', 'v1', 10);
m = AddState(m, 'x3', 'v1', 0);
m = AddState(m, 'x4', 'v1', 0);

%% Add Inputs
m = AddInput(m, 'u1', 'v1', 1);
m = AddInput(m, 'u2', 'v1', @(t,q)2);
m = AddInput(m, 'u3', 'v1', @(t,q)(q), 3);
m = AddInput(m, 'u4', 'v1', @(t,q)(t*q), 4);

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
m = AddReaction(m, 'r1', 'v1', '', '', '', '', 'k1');
m = AddReaction(m, 'r2', '', '', '', 'x3', '', 'k1');
m = AddReaction(m, 'r3', '', '', '', 'x3', 'x4', 'k1');
m = AddReaction(m, 'r4', '', '', '', 'u3', 'x4', 'k1');
m = AddReaction(m, 'r5', '', '', '', 'u3', 'u4', 'k1');
m = AddReaction(m, 'r6', '', '', '', 'u3', '', 'k1');

m = AddReaction(m, 'r7', '', 'x1', '', '', '', 'k1');
m = AddReaction(m, 'r8', '', 'x1', '', 'x3', '', 'k1');
m = AddReaction(m, 'r9', '', 'x1', '', 'x3', 'x4', 'k1');
m = AddReaction(m, 'r10', '', 'x1', '', 'u3', 'x4', 'k1');
m = AddReaction(m, 'r11', '', 'x1', '', 'u3', 'u4', 'k1');
m = AddReaction(m, 'r12', '', 'x1', '', 'u3', '', 'k1');

m = AddReaction(m, 'r13', '', 'x1', 'x2', '', '', 'k1');
m = AddReaction(m, 'r14', '', 'x1', 'x2', 'x3', '', 'k1');
m = AddReaction(m, 'r15', '', 'x1', 'x2', 'x3', 'x4', 'k1');
m = AddReaction(m, 'r16', '', 'x1', 'x2', 'u3', 'x4', 'k1');
m = AddReaction(m, 'r17', '', 'x1', 'x2', 'u3', 'u4', 'k1');
m = AddReaction(m, 'r18', '', 'x1', 'x2', 'u3', '', 'k1');

m = AddReaction(m, 'r19', '', 'x1', 'u2', '', '', 'k1');
m = AddReaction(m, 'r20', '', 'x1', 'u2', 'x3', '', 'k1');
m = AddReaction(m, 'r21', '', 'x1', 'u2', 'x3', 'x4', 'k1');
m = AddReaction(m, 'r22', '', 'x1', 'u2', 'u3', 'x4', 'k1');
m = AddReaction(m, 'r23', '', 'x1', 'u2', 'u3', 'u4', 'k1');
m = AddReaction(m, 'r24', '', 'x1', 'u2', 'u3', '', 'k1');

m = AddReaction(m, 'r25', '', 'u1', 'u2', '', '', 'k1');
m = AddReaction(m, 'r26', '', 'u1', 'u2', 'x3', '', 'k1');
m = AddReaction(m, 'r27', '', 'u1', 'u2', 'x3', 'x4', 'k1');
m = AddReaction(m, 'r28', '', 'u1', 'u2', 'u3', 'x4', 'k1');
m = AddReaction(m, 'r29', '', 'u1', 'u2', 'u3', 'u4', 'k1');
m = AddReaction(m, 'r30', '', 'u1', 'u2', 'u3', '', 'k1');

m = AddReaction(m, 'r31', '', '', 'u2', '', '', 'k1');
m = AddReaction(m, 'r32', '', '', 'u2', 'x3', '', 'k1');
m = AddReaction(m, 'r33', '', '', 'u2', 'x3', 'x4', 'k1');
m = AddReaction(m, 'r34', '', '', 'u2', 'u3', 'x4', 'k1');
m = AddReaction(m, 'r35', '', '', 'u2', 'u3', 'u4', 'k1');
m = AddReaction(m, 'r36', '', '', 'u2', 'u3', '', 'k1');

%% Finalize Model
m = FinalizeModel(m);
