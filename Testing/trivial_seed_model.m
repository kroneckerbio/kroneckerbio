function [m, con, obj, opts] = trivial_seed_model()
m = InitializeModelMassActionAmount('Trivial');
m = AddCompartment(m, 'v', 3, 1);
m = AddState(m, 'x', 'v', 's');
m = AddSeed(m, 's', 1);
m = addStatesAsOutputs(m);
m = FinalizeModel(m);

con = experimentInitialValue(m);

sd = sdLinear(1, 0);
obs = observationLinearWeightedSumOfSquares(1, 1, sd, 'SimpleData');
obj = obs.Objective(2);

opts.Verbose = false;
opts.UseParams = [];
opts.UseSeeds = 1;
end
