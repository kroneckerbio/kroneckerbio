function tests = UT04_Simulation()
tests = functiontests(localfunctions);
end

function testDoseApplication(a)
[m, con, unused, opts] = dose_model();
obs = observationSelect(6);

sim0 = SimulateSystem(m, con(1), obs, opts);

a.verifyEqual(sim0.y(:,end), [0;8;0])

sim1 = SimulateSystem(m, con(2), obs, opts);

a.verifyEqual(sim1.y(:,end), [0;4;0])
end

function testSimulateSimple(a)
[m, con, unused, opts] = simple_model();

sim = SimulateSystem(m, con, 6, opts);

a.verifyEqual(numel(sim.x(4)), m.nx)
a.verifyEqual(numel(sim.u(4)), m.nu)
a.verifyEqual(numel(sim.y(4)), m.ny)
end

function testSimulateSelectSimple(a)
[m, con, unused, opts] = simple_model();

obs = observationSelect([2,4]);

sim = SimulateSystem(m, con, obs, opts);

a.verifyEqual(size(sim.x), [m.nx,2])
a.verifyEqual(size(sim.u), [m.nu,2])
a.verifyEqual(size(sim.y), [m.ny,2])
end

function testSimulateEvent(a)
[m, con, unused, opts] = simple_model();
eve1 = eventDropsBelow(m, 10, 15);
eve2 = eventDropsBelow(m, 1, 2);
obs = observationEvents(6, [eve1;eve2]);

sim = SimulateSystem(m, con, obs, opts);

a.verifyEqual(size(sim.ue,1), m.nu)
a.verifyEqual(size(sim.xe,1), m.nx)
a.verifyEqual(size(sim.ye,1), m.ny)
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

sim = SimulateSystem(m, con, 10, opts);
y = sim.y([2,7]);
a.verifyEqual(size(y), [m.ny,2])
end
