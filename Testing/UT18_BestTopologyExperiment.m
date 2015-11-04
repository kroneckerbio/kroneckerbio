function tests = UT18_BestTopologyExperiment()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testMAPK(a)
[m, con, ~, obj, opts, objPrior] = mapk_models();

target = @entropy;
opts.TargetTol = 0.5; % Bad tolerance to make this go faster
best = BestTopologyExperiment(m, [], [], objPrior, [], [], [], con, obj, target, opts);
end
