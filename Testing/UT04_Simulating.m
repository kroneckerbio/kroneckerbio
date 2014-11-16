function tests = UT04_Simulating()
tests = functiontests(localfunctions);
end

function setupOnce(a)
m = LoadModelMassAction('Simple.txt');
m = addStatesAsOutputs(m);
m = FinalizeModel(m);

a.TestData.m = m;
end

function testSimulateSimple1(a)
m = a.TestData.m;
con = Experiment(m, 6, m.s, false, false, m.u);
sim = Simulate(m,con);
end

function testSimulateSimpleWithDosing(a)
dose = [2;0;0];
schedule = 1:6;
m = a.TestData.m;
con = Experiment(m, 6, m.s, false, false, m.u, dose, schedule);

a.verifyEqual(con.d(0), zeros(m.ns,1))
a.verifyEqual(con.d(1), dose)

sim = Simulate(m,con);
end