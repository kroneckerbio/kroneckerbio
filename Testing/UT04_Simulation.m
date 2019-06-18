function tests = UT04_Simulation()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testDoseApplication(a)
[m, con, ~, opts] = dose_model();
obs = observationSelect(6);

sim0 = SimulateSystem(m, con(1), obs, opts);

a.verifyEqual(sim0.y(:,end), [0;8;0])

sim1 = SimulateSystem(m, con(2), obs, opts);

a.verifyEqual(sim1.y(:,end), [0;4;0])
end

function testSimulateSimple(a)
[m, con, ~, opts] = simple_model();

sim = SimulateSystem(m, con, 6, opts);

a.verifyEqual(numel(sim.x(4)), m.nx)
a.verifyEqual(numel(sim.u(4)), m.nu)
a.verifyEqual(numel(sim.y(4)), m.ny)
end

function testSimulateSimpleAnalytic(a)
[m, con, ~, opts] = simple_analytic_model();

sim = SimulateSystem(m, con, 40, opts);

a.verifyEqual(numel(sim.x(4)), m.nx)
a.verifyEqual(numel(sim.u(4)), m.nu)
a.verifyEqual(numel(sim.y(4)), m.ny)
end

function testSimulateSelectSimple(a)
[m, con, ~, opts] = simple_model();

obs = observationSelect([2,4]);

sim = SimulateSystem(m, con, obs, opts);

a.verifyEqual(size(sim.x), [m.nx,2])
a.verifyEqual(size(sim.u), [m.nu,2])
a.verifyEqual(size(sim.y), [m.ny,2])
end

function testSimulateSimpleSteadyState(a)
simpleopts.steadyState = true;
[m, con, unused, opts, con_sscheck] = simple_model(simpleopts);

sim = SimulateSystem(m, con, 6, opts);

tbig = 10000;
sim_sscheck = SimulateSystem(m, con_sscheck, tbig, opts);

a.verifyEqual(numel(sim.x(4)), m.nx)
a.verifyEqual(numel(sim.u(4)), m.nu)
a.verifyEqual(numel(sim.y(4)), m.ny)

% Verify steady state works correctly by simulating an initial value
% experiment for a long time.
% Dosing adds to states at t=0, which won't be reflected in the steady
% state. Subtract off the dose.
dose_t0 = m.x0(con.d(0)) - m.x0(zeros(m.ns,1));
x_ss_predose = sim.x(0) - dose_t0;
a.verifyEqual(x_ss_predose, sim_sscheck.x(tbig), 'AbsTol', 1e-7)
end

function testSimulateSimpleSteadyStateStartingAtSteadyState(a)

simpleopts.steadyState = true;
[m, con, unused, opts, con_sscheck] = simple_model(simpleopts);

% Set all parameters to zero so that model is inactive and all sets of
% states are steady states
m = m.Update(zeros(size(m.k)));

% Simulate steady state experiment
sim = SimulateSystem(m, con, 6, opts);

% Verify steady state works correctly by comparing to initial value
% Dosing adds to states at t=0, which won't be reflected in the steady
% state. Subtract off the dose.
dose_t0 = m.x0(con.d(0)) - m.x0(zeros(m.ns,1));
x_ss_predose = sim.x(0) - dose_t0;
a.verifyEqual(x_ss_predose, m.x0(m.s), 'AbsTol', 1e-7)
end

function testSimulateEvent(a)
[m, con, ~, opts] = simple_model();
eve1 = eventDropsBelow(m, 10, 10);
eve2 = eventDropsBelow(m, 1, 2);
obs = observationEvents(6, [eve1;eve2]);

sim = SimulateSystem(m, con, obs, opts);

a.verifyEqual(size(sim.ue,1), m.nu)
a.verifyEqual(size(sim.xe,1), m.nx)
a.verifyEqual(size(sim.ye,1), m.ny)
end

function testSimulateTerminalEvent(a)
[m, con, ~, opts] = simple_model();
eve1 = eventDropsBelow(m, 10, 10);
eve2 = eventDropsBelow(m, 1, 3);
obs = observationEvents(6, [eve1;eve2], 1);

sim = SimulateSystem(m, con, obs, opts);

a.verifyNotEqual(sim.te(end), 6);
a.verifyEqual(size(sim.te), [1, 2]);
a.verifyEqual(size(sim.ue), [m.nu, 2])
a.verifyEqual(size(sim.xe), [m.nx, 2])
a.verifyEqual(size(sim.ye), [m.ny, 2])
end

function testSimulateTerminalEventAndAll(a)
[m, con, ~, opts] = simple_model();
eve1 = eventDropsBelow(m, 10, 10);
eve2 = eventDropsBelow(m, 1, 3);
obs1 = observationEvents(6, [eve1;eve2], 1);
obs2 = observationAll(6);

sim = SimulateSystem(m, con, [obs1; obs2], opts);

sim1 = sim(1);

a.verifyNotEqual(sim1.te(end), 6);
a.verifyEqual(size(sim1.te), [1, 2]);
a.verifyEqual(size(sim1.ue), [m.nu, 2])
a.verifyEqual(size(sim1.xe), [m.nx, 2])
a.verifyEqual(size(sim1.ye), [m.ny, 2])

sim2 = sim(2);

a.verifyEqual(sim2.t(end), 6);
end

function testSimulateDifferent(a)
[m, con, ~, opts] = simple_model();
eve1 = eventDropsBelow(m, 10, 15);
eve2 = eventDropsBelow(m, 1, 2);
obs1 = observationEvents(6, [eve1;eve2]);
obs2 = observationSelect(2:6);
obs3 = observationAll(5);
obs = [obs1;obs2;obs3];

sim = SimulateSystem(m, con, obs, opts);

a.verifyEqual(size(sim(1).ue,1), m.nu)
a.verifyEqual(size(sim(1).xe,1), m.nx)
a.verifyEqual(size(sim(1).ye,1), m.ny)

a.verifyEqual(size(sim(2).x), [m.nx,5])
a.verifyEqual(size(sim(2).u), [m.nu,5])
a.verifyEqual(size(sim(2).y), [m.ny,5])

a.verifyEqual(numel(sim(3).x(4)), m.nx)
a.verifyEqual(numel(sim(3).u(4)), m.nu)
a.verifyEqual(numel(sim(3).y(4)), m.ny)
end

function testLinearNoiseApproximationOnSimple(a)
[m, con, ~, opts] = simple_model();
opts.AbsTol = nan;

sim = SimulateLna(m, con, 10, opts);

Vy = sim.Vy([2;4]);
a.verifyEqual(size(Vy), [m.ny^2,2])
end

function testSimulateMichaelisMenten(a)
[m, con, ~, opts] = michaelis_menten_model();

sim = SimulateSystem(m, con, 10, opts);
y = sim.y([2,7]);
a.verifyEqual(size(y), [m.ny,2])
end

function testSimulateEventMichaelisMenten(a)
[m, con, ~, opts, eve] = michaelis_menten_model();

obs = observationEvents(10, eve);

sim = SimulateSystem(m, con, obs, opts);

a.verifyEqual(size(sim.ue,1), m.nu)
a.verifyEqual(size(sim.xe,1), m.nx)
a.verifyEqual(size(sim.ye,1), m.ny)
end

function testDoseList(a)
[m, ~, ~, opts] = dose_model();

dos = doseList(m, [4;3;2], [2;4;5], [2;1;1]);
con = experimentInitialValue(m, [], [], dos);

obs = observationSelect([4, 6]);

sim0 = SimulateSystem(m, con, obs, opts);

a.verifyEqual(sim0.y, [0, 0; 3, 5; 4 4])
end

function testIntegrationFailure(a)
m = InitializeModelAnalytic('infinity');
m = AddCompartment(m, 'v', 3, 1);
m = AddState(m, 'x', 'v', 4);
m = AddReaction(m, '', 'x', [], '1.5'); % Malicious reaction
m = FinalizeModel(m);

con = experimentInitialValue(m);

opts.Verbose = false;

state = warning;
finished = onCleanup(@() warning(state));
warning('off', 'MATLAB:ode15s:IntegrationTolNotMet')
a.verifyError(@()SimulateSystem(m, con, 10, opts), 'KroneckerBio:accumulateOde:IntegrationFailure');
end

function testSymbolicInputArgumentsMassAction(a)
m = simple_model();

t = sym('t');
x = sym('x', [m.nx,1]);
u = sym('u', [m.nu,1]);
k = sym('k', [m.nk,1]);
msym = m.Update(k);

tval = rand;
xval = rand(size(x));
uval = rand(size(u));
kval = rand(size(k));
mval = m.Update(kval);

funs = {'f', 'dfdx', 'dfdu', 'dfdk', 'r', 'drdx', 'drdu', 'drdk', 'v', ...
    'dvdx', 'dvdu', 'y', 'dydx', 'dydu', 'dydk'};
for i = 1:numel(funs)
    val_sym = msym.(funs{i})(t, x, u);
    val = double(subs(val_sym, [t;x;u;k], [tval;xval;uval;kval]));
    val_test = mval.(funs{i})(tval, xval, uval);
    abs_err = abs(val - val_test);
    if ~isempty(abs_err)
        a.verifyLessThan(abs_err, 1e-11)
    end
end
end

function testSymbolicInputArgumentsAnalytic(a)
m = simple_analytic_model();

t = sym('t');
x = sym('x', [m.nx,1]);
u = sym('u', [m.nu,1]);
k = sym('k', [m.nk,1]);
msym = m.Update(k);

tval = rand;
xval = rand(size(x));
uval = rand(size(u));
kval = rand(size(k));
mval = m.Update(kval);

funs = {'f', 'dfdx', 'dfdu', 'dfdk', 'r', 'drdx', 'drdu', 'drdk', ...
    'y', 'dydx', 'dydu', 'dydk'};
for i = 1:numel(funs)
    val_sym = msym.(funs{i})(t, x, u);
    val = double(subs(val_sym, [t;x;u;k], [tval;xval;uval;kval]));
    val_test = mval.(funs{i})(tval, xval, uval);
    abs_err = abs(val - val_test);
    if ~isempty(abs_err)
        a.verifyLessThan(abs_err, 1e-11)
    end
end
end
