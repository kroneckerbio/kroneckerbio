function tests = UT10_Fitting()
tests = functiontests(localfunctions);
end

function setupOnce(a)
m = LoadModelMassAction('Simple.txt');
m = addStatesAsOutputs(m);
m = FinalizeModel(m);

schedule = 1:6;
[d, dddq] = dosingConstant(m.ns, 1, schedule, 1);
con = Experiment(m, 6, m.s, false, false, m.u, d, schedule, 2, [], dddq);

obj = generateFakeData(m, con, 1:m.ny, linspace(0,6,11), sdLinear(0.01,0.1), true, true, 1337);

a.TestData.m = m;
a.TestData.con = con;
a.TestData.obj = obj;
end

function testSimpleFitting(a)
m = a.TestData.m;
con = a.TestData.con;
obj = a.TestData.obj;

opts.UseParams = 1:m.nk;
opts.UseSeeds = 1:m.ns;
opts.UseControls = 1:con.nq;
opts.MaxIter = 2;

Gold = ObjectiveValue(m, con, obj, opts);

[m, con] = FitObjective(m, con, obj, opts);

Gnew = ObjectiveValue(m, con, obj, opts);

a.verifyLessThan(Gnew, Gold)
end