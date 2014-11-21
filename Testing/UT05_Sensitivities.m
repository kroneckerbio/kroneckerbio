function tests = UT05_Sensitivities()
tests = functiontests(localfunctions);
end

function setupOnce(a)
m = LoadModelMassAction('Simple.txt');
m = addStatesAsOutputs(m);
m = FinalizeModel(m);

a.TestData.m = m;
end

function con = dosingExperiment(a)
m = a.TestData.m;
schedule = 1:6;
[d, dddq] = dosingConstant(m.ns, 1, schedule, 1);
con = Experiment(m, 6, m.s, false, false, m.u, d, schedule, 2, [], dddq);
end

function testDosingConstant(a)
m = a.TestData.m;
con = dosingExperiment(a);

a.verifyEqual(con.d(0), zeros(m.ns,1))
a.verifyEqual(con.d(1), [2;0;0])

a.verifyEqual(con.dddq(0), zeros(m.ns,1))
a.verifyEqual(con.dddq(1), [1;0;0])
end

function testSimulateSimpleWithDosing(a)
m = a.TestData.m;
con = dosingExperiment(a);

opts.UseParams = [];
opts.UseSeeds = [];
opts.UseControls = 1;

sim = SimulateSensitivity(m,con,opts);

dydT = sim.dydT(6);
a.verifyEqual(numel(dydT), m.ny)
end

function testSimulateSimpleWithDosingAndParameters(a)
m = a.TestData.m;
con = dosingExperiment(a);

opts.UseParams = 1:m.nk;
opts.UseSeeds = [];
opts.UseControls = 1;

sim = SimulateSensitivity(m,con,opts);

dydT = sim.dydT(6);
a.verifyEqual(numel(dydT), m.ny*(1+m.nk))
end

function testSimulateSelectSimpleWithDosingAndParameters(a)
m = a.TestData.m;
con = dosingExperiment(a);

opts.UseParams = 1:m.nk;
opts.UseSeeds = [];
opts.UseControls = 1;

sim = SimulateSensitivitySelect(m,con,linspace(0,6,101),opts);

dydT = sim.dydT;
a.verifyEqual(size(dydT), [m.ny*(1+m.nk), 101])
end
