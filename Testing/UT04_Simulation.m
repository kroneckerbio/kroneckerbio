function tests = UT04_Simulation()
tests = functiontests(localfunctions);
end

function testDoseApplication(a)
[m, con, unused, opts] = dose_model();

sim0 = SimulateSelect(m, con(1), 6, opts);

a.verifyEqual(sim0.y(:,end), [0;8;0])

sim1 = SimulateSelect(m, con(2), 6, opts);

a.verifyEqual(sim1.y(:,end), [0;4;0])
end

function testSimulateSimple(a)
[m, con, unused, opts] = simple_model();

sim = Simulate(m, con, opts);

y = sim.y(4);
a.verifyEqual(numel(y), m.ny)
end

function testSimulateSelectSimple(a)
[m, con, unused, opts] = simple_model();

sim = SimulateSelect(m, con, [2,4], opts);

y = sim.y;
a.verifyEqual(size(y), [m.ny,2])
end

function testLinearNoiseApproximationOnSimple(a)
[m, con, unused, opts] = simple_model();
opts.AbsTol = nan;

sim = SimulateLna(m, con, opts);

Vy = sim.Vy([2;4]);
a.verifyEqual(size(Vy), [m.ny^2,2])
end

function testMichaelisMenten(a)
[m, con, unused, opts] = michaelis_menten_model();

sim = Simulate(m,con,opts);
y = sim.y([2,7]);
a.verifyEqual(size(y), [m.ny,2])
end
