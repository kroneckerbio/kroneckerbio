function tests = UT08_Gradients()
tests = functiontests(localfunctions);
end

function setupOnce(a)
m = LoadModelMassAction('Simple.txt');
m = addStatesAsOutputs(m);
m = FinalizeModel(m);

schedule = 1:6;
[d, dddq] = dosingConstant(m.ns, 1, schedule, 1);
con = Experiment(m, 6, m.s, false, false, m.u, d, schedule, 2, [], dddq);

obj = generateFakeData(m, con, 1:m.ny, linspace(0,6,11), sdLinear(0.01,0.1), false, true, 42);

a.TestData.m = m;
a.TestData.con = con;
a.TestData.obj = obj;
end

function testObjectiveGradientForward(a)
m = a.TestData.m;
con = a.TestData.con;
obj = a.TestData.obj;

opts.UseParams = 1:m.nk;
opts.UseSeeds = 1:m.ns;
opts.UseControls = 1:con.nq;

% Forward
opts.UseAdjoint = false;
Dfwd = ObjectiveGradient(m, con, obj, opts);

a.verifyEqual(size(Dfwd), [m.nk+m.ns+con.nq,1])

% Adjoint
opts.UseAdjoint = true;
Dadj = ObjectiveGradient(m, con, obj, opts);

a.verifyEqual(size(Dadj), [m.nk+m.ns+con.nq,1])

% Discrete
Ddisc = FiniteObjectiveGradient(m, con, obj, opts);

a.verifyEqual(size(Ddisc), [m.nk+m.ns+1,1])

% THE MOMENT OF TRUTH
a.verifyEqual(Dfwd, Ddisc, 'RelTol', 0.01)
a.verifyEqual(Dadj, Ddisc, 'RelTol', 0.01)
end
