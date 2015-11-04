function tests = UT16_TopologyProbability()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testMAPK(a)
[m, con, ~, obj, opts, objPrior] = mapk_models();

pmy = TopologyProbability(m, con, obj, objPrior, [], [], [], opts);
a.verifyGreaterThan(pmy(1), 0.05)
end
