%% T02 Modifying Models
% After model finalization, components can be removed and additional components
%   can be added.
% Remember to finalized the modified model.

%% Build initial model
m = InitializeModelMassActionAmount('Testing');

m = AddCompartment(m, 'v1', 3, 1);
m = AddCompartment(m, 'v2', 3, 1);

m = AddState(m, 'x1', 'v1', 1);
m = AddState(m, 'x2', 'v1', 0);
m = AddState(m, 'x2', 'v2', 0);

m = AddInput(m, 'u1', 'v1', 1);

m = AddOutput(m, 'y1', 'v1.x2');
m = AddOutput(m, 'y2', {'v1.x2', 'v2.x2'});

m = AddParameter(m, 'k1', 1);
m = AddParameter(m, 'k2', 2);

m = AddReaction(m, 'r1', {'u1', 'x1'}, {'x2'}, 'k1', '', 'v1');
m = AddReaction(m, 'r2', {'v1.x2'}, {'v2.x2'}, 'k2');

% Finialize initial model
m = FinalizeModel(m);

%% Modify model
m = RemoveCompartment(m, 'v2');

m = RemoveState(m, 'v2.x2');

m = RemoveParameter(m, 'k2');

m = RemoveReaction(m, 'r2');

m = RemoveOutput(m, 'y2');

% Finialize modified model
m = FinalizeModel(m);