function [m, con, obj, opts] = equilibrium_model()
m = LoadModel('Equilibrium.txt');

con = experimentInitialValue(m);

sd = sdLinear(0, 1);
values = [
    3, 1, 0.3;
    ];
obs = observationLogWeightedSumOfSquares(values(:,1), values(:,2), sd);
obj = obs.Objective(values(:,3));

opts.Verbose = false;
opts.UseParams = 1;
opts.UseSeeds = [];
opts.UseInputControls = [];
opts.UseDoseControls = [];
end
