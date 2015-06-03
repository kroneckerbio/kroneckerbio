function tests = UT11_Probability()
tests = functiontests(localfunctions);
end

function testObjectiveProbabilitySimple(a)
[m, con, obj, opts] = simple_model();

p = ObjectiveProbability(m, con, obj, opts);
end

function testObjectiveLogLikelihoodSimple(a)
[m, con, obj, opts] = simple_model();

logp = ObjectiveLogLikelihood(m, con, obj, opts);
end

function testObjectiveProbabilityMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();

p = ObjectiveProbability(m, con, obj, opts);
end

function testObjectiveLogLikelihoodMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();

logp = ObjectiveLogLikelihood(m, con, obj, opts);
end