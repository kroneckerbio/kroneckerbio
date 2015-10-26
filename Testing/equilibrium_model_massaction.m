function m = equilibrium_model_massaction
% Basic massaction model with seeds
m = MassActionAmountModel('Equilibrium_MassAction');

m.AddCompartment('Solution', 3, 1);

m.AddSeed('A_0', 1);
m.AddSeed('B_0', 2);
m.AddSeed('C_0', 0);

m.AddState('A', 'Solution', 'A_0');
m.AddState('B', 'Solution', 'B_0');
m.AddState('C', 'Solution', 'C_0');

m.AddOutput('A', 'A');
m.AddOutput('B', 'B');
m.AddOutput('C', 'C');

m.AddParameter('kf', 5);
m.AddParameter('kr', 3);

m.AddReaction('r1', {'A', 'B'}, {'C'}, 'kf', 'kr');

m.Finalize;