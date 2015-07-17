function tests = UT08_Gradient()
tests = functiontests(localfunctions);
end

function testDoseModel(a)
[m, con, obj, opts] = dose_model();

% Forward
opts.UseAdjoint = false;
Dfwd = ObjectiveGradient(m, con(1), obj, opts);

% Adjoint
opts.UseAdjoint = true;
Dadj = ObjectiveGradient(m, con(1), obj, opts);

% Discrete
Ddisc = FiniteObjectiveGradient(m, con(1), obj, opts);

a.verifyEqual(Dfwd, Ddisc, 'RelTol', 0.001)
a.verifyEqual(Dfwd, Dadj, 'RelTol', 0.001)

% Forward
opts.UseAdjoint = false;
Dfwd = ObjectiveGradient(m, con(2), obj, opts);

% Adjoint
opts.UseAdjoint = true;
Dadj = ObjectiveGradient(m, con(2), obj, opts);

% Discrete
Ddisc = FiniteObjectiveGradient(m, con(2), obj, opts);

a.verifyEqual(Dfwd, Ddisc, 'RelTol', 0.001)
a.verifyEqual(Dfwd, Dadj, 'RelTol', 0.001)

% Forward
opts.UseAdjoint = false;
Dfwd = ObjectiveGradient(m, con(3), obj, opts);

% Adjoint
opts.UseAdjoint = true;
Dadj = ObjectiveGradient(m, con(3), obj, opts);

% Discrete
Ddisc = FiniteObjectiveGradient(m, con(3), obj, opts);

a.verifyEqual(Dfwd, Ddisc, 'RelTol', 0.001)
a.verifyEqual(Dfwd, Dadj, 'RelTol', 0.001)
end

function testObjectiveGradientSimpleSteadyState(a)
simpleopts.steadyState = true;
[m, con, obj, opts] = simple_model(simpleopts);
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

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

a.verifyEqual(Dfwd, Ddisc, 'RelTol', 0.001, 'AbsTol', 0.001)
a.verifyEqual(Dadj, Ddisc, 'RelTol', 0.001, 'AbsTol', 0.001)
end

function testObjectiveGradientSimple(a)
[m, con, obj, opts] = simple_model();
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

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

a.verifyEqual(Dfwd, Ddisc, 'RelTol', 0.001, 'AbsTol', 0.001)
a.verifyEqual(Dadj, Ddisc, 'RelTol', 0.001, 'AbsTol', 0.001)
end

function testObjectiveGradientMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

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

a.verifyEqual(Dfwd, Ddisc, 'RelTol', 0.001)
a.verifyEqual(Dadj, Ddisc, 'RelTol', 0.001)
end