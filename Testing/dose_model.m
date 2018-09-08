function [m, con, obj, opts] = dose_model()
m = LoadModel('DoseModel.txt');
m = addStatesAsOutputs(m);
m = FinalizeModel(m);

dos0 = doseConstant(m, 2, 0:3, 1);
con0 = experimentInitialValue(m, [], [], dos0);

dos1 = doseConstant(m, 2, 5:7, 1);
con1 = experimentInitialValue(m, [], [], dos1);

dos2 = doseConstant(m, 3, [3,5], 2);
con2 = experimentInitialValue(m, [], 0.5, dos2);

con = [con0;con1;con2];

sd = sdLinear(0.1, 1);
values = [
    1, 4, 2;
    2, 5, 5;
    3, 5.5, 3;
    ];
obs = observationLinearWeightedSumOfSquares(values(:,1), values(:,2), sd, 'SimpleData');
obj = obs.Objective(values(:,3));

opts.Verbose = false;
opts.UseParams = 1;
opts.UseSeeds = [];
opts.UseInputControls = [];
opts.UseDoseControls = [];
end
