%% T01a Building a Massaction Model
% Build a massaction model, demonstrating different ways of calling the
%   component-adding functions.
% Components can be added in any order. Checks are done when the model is
%   finalized.
% For detailed usage, see the documentation in the source files for each
%   function.

%% Initialize Model
m = InitializeModelMassActionAmount('Testing');

%% Add Compartments
% Each species (state and input) must belong to a compartment
% All different dimensions
m = AddCompartment(m, 'v1', 3, 1);
m = AddCompartment(m, 'v2', 2, 1);
m = AddCompartment(m, 'v3', 1, 1);
m = AddCompartment(m, 'v4', 0, 1);

%% Add Seeds
% Seeds are special parameters that control state initial amounts. State initial
%   amounts can be linear combinations of seeds.
m = AddSeed(m, 's1', 5);

%% Add States
% States are the usual species that react and change over time, the terms that
%   appear as time derivatives.
% The 4th argument can be blank (default 0 initial amount), a number for initial
%   amount, or a linear combination of seeds (see AddStateMassActionAmount documentation)
% State names must be unique per compartment but a state can appear in multiple
%   compartments. In reaction and output specifications, the state name can be
%   given as the unqualified 'state' or the qualified 'compartment.state'.
m = AddState(m, 'x1', 'v1');
m = AddState(m, 'x2', 'v1', 10);
m = AddState(m, 'x3', 'v1', 0);
m = AddState(m, 'x4', 'v1', 's1');

%% Add Inputs
% Inputs are mathematical functions that change with time. A common use is for
%   species that don't change amount during the simulation, such as a substrate
%   in large excess.
% A constant baseline input amount is specified here. The time-varying or more
%   complicated behavior is added later in experiments.
% Note that inputs are specified somewhat like states here and require a
%   compartment.
% Input names must be unique per compartment but a state can appear in multiple
%   compartments. In reaction and output specifications, the input name can be
%   given as the unqualified 'input' or the qualified 'compartment.input'.
m = AddInput(m, 'u1', 'v1', 1);
m = AddInput(m, 'u2', 'v1', 2);
m = AddInput(m, 'u3', 'v1', 3);
m = AddInput(m, 'u4', 'v1', 4);

%% Add Outputs
% Outputs are quantities of interest to the user, either in plotting their
%   simulated amounts over time, or in comparing to data during fitting.
% In massaction models, outputs are linear combinations of states and inputs.
%   See the addOutputMassActionAmount documentation.
m = AddOutput(m, 'y1', 'x1');
m = AddOutput(m, 'y2', {'x1', 2});
m = AddOutput(m, 'y3', {'u1', 3});
m = AddOutput(m, 'y4', {'x1', 1; 'u1', 2});
m = addOutputAsRegex(m, 'y5', {'x'}); % optional method of specifying an output as a linear combo of regex matches
m = addOutputAsRegex(m, 'y6', {'u', 2});

%% Add Parameters
% Parameters are the rate constants in massaction models
m = AddParameter(m, 'k1', 1);
m = AddParameter(m, 'k2', 2);
m = AddParameter(m, 'k3', 3);

%% Add Reactions
% Reactions change the amounts of species over time. In massaction models,
%   reactions can have 0, 1, or 2 reactants and products.
% Optional arguments that appear in the middle of the arguments list are
%   specified by anything that returns isempty == true, like [], {}, or ''.
% The 2nd argument, the reaction name, is optional for model building but
%   required for later model modification.
% The last argument, the reaction compartment, is optional but specifies the
%   reaction the unqualified species (states and inputs) are in.
% In massaction models, a warning will be thrown if identical reactions (same
%   name, reactants, products, rate constants) are added.
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
m = AddReaction(m, 'r37', 'x1', {'x1', 'x1'}, 'k1', '', 'v1');

% Add identical reaction - throws warning on model finalization
m = AddReaction(m, 'r22', {'x1', 'u2'}, {'x4', 'u3'}, 'k1');

%% Finalize Model
% Preparing the model for use with the rest of kroneckerbio (simulation,
%   fitting, analysis) is relatively expensive, with the calculation of different
%   forms of derivatives. So it's done once after all desired components are
%   added.
% Components can be added or removed after finalization, but re-finalization is
%   necessary.
% Remember to finalize models before use.
m = FinalizeModel(m);
