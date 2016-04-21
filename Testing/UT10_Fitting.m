function tests = UT10_Fitting()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testSimpleFitting(a)
[m, con, obj, opts] = simple_model();
opts.MaxIter = 2;

Gold = ObjectiveValue(m, con, obj, opts);

[m, con] = FitObjective(m, con, obj, opts);

Gnew = ObjectiveValue(m, con, obj, opts);

a.verifyLessThan(Gnew, Gold)
end

function testSimpleAnalyticFitting(a)
[m, con, obj, opts] = simple_analytic_model();
opts.MaxIter = 2;

Gold = ObjectiveValue(m, con, obj, opts);

[m, con] = FitObjective(m, con, obj, opts);

Gnew = ObjectiveValue(m, con, obj, opts);

a.verifyLessThan(Gnew, Gold)
end

function testMichaelisMentenFitting(a)
[m, con, obj, opts] = michaelis_menten_model();

Gold = ObjectiveValue(m, con, obj, opts);

timer = tic;
[m, con, G, D] = FitObjective(m, con, obj, opts);
time = toc(timer);

Gnew = ObjectiveValue(m, con, obj, opts);
Dnew = ObjectiveGradient(m, con, obj, opts);

a.verifyLessThan(Gnew, Gold)
a.verifyLessThan(time, 12)

a.verifyEqual(G, Gnew, 'RelTol', 1e-4)
a.verifyEqual(D, Dnew, 'RelTol', 1e-4)
end
