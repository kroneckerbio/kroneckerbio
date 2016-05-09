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

function testObjectiveHessianMichaelisMentenMultipleExperiments(a)
[m, con, obj, opts] = michaelis_menten_model();

con = repmat(con,2,1);
obj = repmat(obj,1,2);
opts.UseSeeds = repmat(logical(opts.UseSeeds),1,2);
opts.UseSeeds(opts.UseSeeds(:,2),2) = false(size(opts.UseSeeds(:,2)));
opts.UseInputControls = repmat({logical(opts.UseInputControls)},2,1);
opts.UseInputControls{2} = false(size(opts.UseInputControls{1}));
opts.UseDoseControls = repmat({logical(opts.UseDoseControls)},2,1);
opts.UseDoseControls{2} = false(size(opts.UseDoseControls{1}));
opts.AbsTol = 1e-9;

verifyHessian(a, m, con, obj, opts)
end

function verifyHessian(a, m, con, obj, opts)
opts.ComplexStep = true;
nTk = nnz(opts.UseParams);
nTs = nnz(opts.UseSeeds);
if iscell(opts.UseInputControls)
    nTq = sum(cellfun(@nnz, opts.UseInputControls));
else
    nTq = nnz(opts.UseInputControls);
end
if iscell(opts.UseDoseControls)
    nTh = sum(cellfun(@nnz, opts.UseDoseControls));
else
    nTh = nnz(opts.UseDoseControls);
end
nT = nTk + nTs + nTq + nTh;

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
