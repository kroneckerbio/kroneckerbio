function tests = UT17_BestParameterExperiment()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testDoseModel(a)
[m, con, obs, opts] = dose_model();
obs = repmat(obs, 1,numel(con));

goal = @(F)goalParameterSpaceSumOfLog(F);
best = BestParameterExperiment(m, [], [], con, obs, goal, opts);
end
