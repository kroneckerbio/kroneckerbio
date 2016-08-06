function tests = UT08_Gradient()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testDoseModel(a)
[m, con, obj, opts] = dose_model();

verifyGradient(a, m, con(1), obj, opts)
verifyGradient(a, m, con(1), obj, opts)
verifyGradient(a, m, con(1), obj, opts)
end

function testHigherOrderDose(a)
[m, con, obj, opts] = higher_order_dose_model();

verifyGradient(a, m, con, obj, opts)
end

function testObjectiveGradientSimpleSteadyState(a)
simpleopts.steadyState = true;
[m, con, obj, opts] = simple_model(simpleopts);

verifyGradient(a, m, con, obj, opts)
end

function testObjectiveGradientSimpleSteadyStateStartingAtSteadyState(a)
simpleopts.steadyState = true;
[m, con, obj, opts] = simple_model(simpleopts);

% Set k, q, and s to zero so that model is inactive in the basal
% simulation and the starting state is a steady state
newk = repmat(1e-9, size(m.k));
m = m.Update(newk);
newq = repmat(1e-9, size(con.q));
news = repmat(1e-9, size(con.s));
con = con.Update(news, newq, con.h);

verifyGradient(a, m, con, obj, opts)
end

function testObjectiveGradientSimple(a)
[m, con, obj, opts] = simple_model();

verifyGradient(a, m, con, obj, opts)
end

function testObjectiveGradientSimpleKineticPrior(a)
simpleopts.kineticPrior = true;
[m, con, obj, opts] = simple_model(simpleopts);
obj = obj(2);

verifyGradient(a, m, con, obj, opts)
end

function testObjectiveGradientSimpleLogKineticPrior(a)
simpleopts.logKineticPrior = true;
[m, con, obj, opts] = simple_model(simpleopts);
obj = obj(2);

verifyGradient(a, m, con, obj, opts)
end

function testObjectiveGradientSimpleLogSeedPrior(a)
simpleopts.logSeedPrior = true;
[m, con, obj, opts] = simple_model(simpleopts);
obj = obj(2);

verifyGradient(a, m, con, obj, opts)
end

function testObjectiveGradientSimpleAnalytic(a)
[m, con, obj, opts] = simple_analytic_model();

verifyGradient(a, m, con, obj, opts)
end

function testObjectiveGradientMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();

verifyGradient(a, m, con, obj, opts)
end

function verifyGradient(a, m, con, obj, opts)
opts.ComplexStep = true;
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

opts.Normalized = false;
% Forward
opts.UseAdjoint = false;
Dfwd = ObjectiveGradient(m, con, obj, opts);

% Adjoint
opts.UseAdjoint = true;
Dadj = ObjectiveGradient(m, con, obj, opts);

% Discrete
Ddisc = FiniteObjectiveGradient(m, con, obj, opts);

a.verifyEqual(size(Dfwd), [nT,1])
a.verifyEqual(size(Dadj), [nT,1])
a.verifyEqual(size(Ddisc), [nT,1])

a.verifyEqual(Dfwd, Ddisc, 'RelTol', 0.001, 'AbsTol', 1e-4)
a.verifyEqual(Dadj, Ddisc, 'RelTol', 0.001, 'AbsTol', 1e-4)

opts.Normalized = true;
% Forward
opts.UseAdjoint = false;
Dfwdn = ObjectiveGradient(m, con, obj, opts);

% Adjoint
opts.UseAdjoint = true;
Dadjn = ObjectiveGradient(m, con, obj, opts);

% Discrete
Ddiscn = FiniteObjectiveGradient(m, con, obj, opts);

a.verifyEqual(Dfwdn, Ddiscn, 'RelTol', 0.001, 'AbsTol', 1e-4)
a.verifyEqual(Dadjn, Ddiscn, 'RelTol', 0.001, 'AbsTol', 1e-4)
end
