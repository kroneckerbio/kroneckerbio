function tests = UT09_Hessian()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testHigherOrderDose(a)
[m, con, obj, opts] = higher_order_dose_model();

verifyHessian(a, m, con, obj, opts)
end

function testObjectiveHessianSimple(a)
[m, con, obj, opts] = simple_model();

verifyHessian(a, m, con, obj, opts)
end

function testDoseModel(a)
[m, con, obj, opts] = dose_model();

verifyHessian(a, m, con(1), obj, opts)
verifyHessian(a, m, con(2), obj, opts)
verifyHessian(a, m, con(3), obj, opts)
end

function testObjectiveHessianSimpleSteadyState(a)
simpleopts.steadyState = true;
[m, con, obj, opts] = simple_model(simpleopts);

verifyHessian(a, m, con, obj, opts)
end

function testObjectiveHessianSimpleKineticPriors(a)
simpleopts.kineticPrior = true;
[m, con, obj, opts] = simple_model(simpleopts);
obj = obj(2);

verifyHessian(a, m, con, obj, opts)
end

function testObjectiveHessianSimpleLogKineticPriors(a)
simpleopts.logKineticPrior = true;
[m, con, obj, opts] = simple_model(simpleopts);
obj = obj(2);

verifyHessian(a, m, con, obj, opts)
end

function testObjectiveHessianSimpleLogSeedPriors(a)
simpleopts.logSeedPrior = true;
[m, con, obj, opts] = simple_model(simpleopts);
obj = obj(2);

verifyHessian(a, m, con, obj, opts)
end

function testObjectiveHessianMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();

verifyHessian(a, m, con, obj, opts)
end

function testObjectiveHessianSimpleAnalytic(a)
[m, con, obj, opts] = simple_analytic_model();

verifyHessian(a, m, con, obj, opts)
end

function verifyHessian(a, m, con, obj, opts)
opts.ComplexStep = true;
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

opts.Normalized = false;
Hfwd = ObjectiveHessian(m, con, obj, opts);
Hdisc = FiniteObjectiveHessian(m, con, obj, opts);

a.verifyEqual(size(Hfwd), [nT,nT])
a.verifyEqual(size(Hdisc), [nT,nT])

a.verifyEqual(Hfwd, Hdisc, 'RelTol', 0.001, 'AbsTol', 1e-4)

opts.Normalized = true;
Hfwdn = ObjectiveHessian(m, con, obj, opts);
Hdiscn = FiniteObjectiveHessian(m, con, obj, opts);

a.verifyEqual(Hfwdn, Hdiscn, 'RelTol', 0.001, 'AbsTol', 1e-4)
end
