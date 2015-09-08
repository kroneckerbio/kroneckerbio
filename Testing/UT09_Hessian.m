function tests = UT09_Hessian()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testObjectiveHessianSimple(a)
[m, con, obj, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

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
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

verifyHessian(a, m, con, obj, opts)
end

function testObjectiveHessianMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

verifyHessian(a, m, con, obj, opts)
end

function testObjectiveHessianSimpleAnalytic(a)
[m, con, obj, opts] = simple_analytic_model();
%m = m.Update((rand(m.nk,1)+0.5)./2);
%con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

verifyHessian(a, m, con, obj, opts)
end

function verifyHessian(a, m, con, obj, opts)
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

Hfwd = ObjectiveHessian(m, con, obj, opts);

opts.ImaginaryStep = true;
Hdisc = FiniteObjectiveHessian(m, con, obj, opts);

a.verifyEqual(size(Hfwd), [nT,nT])
a.verifyEqual(size(Hdisc), [nT,nT])

a.verifyEqual(Hfwd, Hdisc, 'RelTol', 0.001, 'AbsTol', 1e-4)
end