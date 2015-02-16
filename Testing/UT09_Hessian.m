function tests = UT09_Hessian()
tests = functiontests(localfunctions);
end

function testObjectiveGradientSimple(a)
[m, con, obj, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

% Forward
Hfwd = ObjectiveHessian(m, con, obj, opts);

% Discrete
Hdisc = FiniteObjectiveHessian(m, con, obj, opts);

a.verifyEqual(size(Hfwd), [nT,nT])
a.verifyEqual(size(Hdisc), [nT,nT])

a.verifyEqual(Hfwd, Hdisc, 'RelTol', 0.001)
a.verifyEqual(Dadj, Hdisc, 'RelTol', 0.001)
end
