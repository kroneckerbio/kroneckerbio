function tests = UT07_Objective()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testObjectiveValueSimple(a)
[m, con, obj, opts] = simple_model();

G = ObjectiveValue(m, con, obj, opts);
end

function testHigherOrderDose(a)
[m, con, obj, opts] = higher_order_dose_model();

obs = observationSelect(0:10);
sim = SimulateSystem(m, con, obs, opts);

int = sim.int;
t = 0;
verifyDerivatives(a, obj, int, t)
end

function testObjectiveWeightedSumOfSquares(a)
[m, con, obj, opts] = simple_model(); % objectiveWeightedSumOfSquares is the default obj fun

obs = observationSelect(1:6);
sim = SimulateSystem(m, con, obs, opts);

int = sim.int;
t = 1;
verifyDerivatives(a, obj, int, t)
end

function testObjectiveWeightedSumOfSquaresSteadyState(a)
simpleopts.steadyState = true;
[m, con, obj, opts] = simple_model(simpleopts); % objectiveWeightedSumOfSquares is the default obj fun

obs = observationSelect(1:6);
sim = SimulateSystem(m, con, obs, opts);

int = sim.int;
t = 1;
verifyDerivatives(a, obj, int, t)
end

function testObjectiveValueSimpleAnalytic(a)
[m, con, obj, opts] = simple_analytic_model();

G = ObjectiveValue(m, con, obj, opts);
end

function testObjectiveWeightSumOfSquaresDose(a)
[m, con, obj, opts] = dose_model();

sim = SimulateSystem(m, con(1), obj, opts);

int = sim.int;
t = 4;
verifyDerivatives(a, obj, int, t)
end

function testObjectiveValueMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();

G = ObjectiveValue(m, con, obj, opts);
end

function testObjectiveWeightedSumOfSquaresMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();

obs = observationSelect(1:10);
sim = SimulateSystem(m, con, obs, opts);

int = sim.int;
t = 1;
verifyDerivatives(a, obj, int, t);

end

function testObjectiveLogWeightedSumOfSquares(a)
simpleopts.objectiveFun = 'observationLogWeightedSumOfSquares';
[m, con, obj, opts] = simple_model(simpleopts);

obs = observationSelect(1:6);
sim = SimulateSystem(m, con, obs, opts);

int = sim.int;
t = 1;
verifyDerivatives(a, obj, int, t)
end

function testObjectiveLogWeightSumOfSquaresDose(a)
[m, con, ~, opts] = dose_model();

obs = observationLogWeightedSumOfSquares(2, 5, sdLinear(0.2, 2));
obj = obs.Objective(4);
sim = SimulateSystem(m, con(1), obs, opts);

int = sim.int;
t = 5;
verifyDerivatives(a, obj, int, t)
end

% function testObjectiveWeightedSumOfSquaresNonNeg(a)
% simpleopts.objectiveFun = 'objectiveWeightedSumOfSquaresNonNeg';
% [m, con, obj, opts] = simple_model(simpleopts);
% 
% obs = observationSelect(1:6);
% sim = SimulateSystem(m, con, obs, opts);
% 
% int = sim.int;
% t = 1;
% verifyDerivatives(a, obj, int, t)
% end

function testObjectiveNormalPriorOnKineticParameters(a)
[m, con, ~, opts] = simple_model(); % objectiveWeightedSumOfSquares is the default obj fun

nk = m.nk;
nTk = sum(opts.UseParams);
kbar = rand(nk,1);
Vbar = rand(nk);
Vbar = Vbar*Vbar.'; % Ensure symmetric and positive definite

obj = objectiveNormalPriorOnKineticParameters(kbar, Vbar, 'Testing kinetic priors');

% Check objective value
G = ObjectiveValue(m, con, obj, opts);
deltaTk = m.k-kbar; 
deltaTk = deltaTk(opts.UseParams);
G_expected = deltaTk.'/Vbar(opts.UseParams,opts.UseParams)*deltaTk;
a.verifyEqual(G, G_expected, 'AbsTol', 1e-9, 'RelTol', 1e-6)

% Check derivatives wrt parameters using finite differences
sim = SimulateSystem(m, con, 0, opts);
int = sim.int;
int.UseParams = opts.UseParams; % Add UseParams field
t = 0;
verifyDerivativesParameters(a, obj, int)
end

function testObjectiveLogNormalPriorOnKineticParameters(a)
[m, con, ~, opts] = simple_model(); % objectiveWeightedSumOfSquares is the default obj fun

nk = m.nk;
nTk = sum(opts.UseParams);
kbar = rand(nk,1);
Vlogbar = rand(nk);
Vlogbar = Vlogbar*Vlogbar.'; % Ensure symmetric and positive definite
Vbar = diag(kbar)*Vlogbar*diag(kbar); % Multiply by square parameters because objective function takes "non-normalized" inputs

obj = objectiveLogNormalPriorOnKineticParameters(kbar, Vbar, 'Testing log kinetic priors');

% Check objective value
G = ObjectiveValue(m, con, obj, opts);
deltaTk = log(m.k)-log(kbar); 
deltaTk = deltaTk(opts.UseParams);
G_expected = deltaTk.'/Vlogbar(opts.UseParams,opts.UseParams)*deltaTk;
a.verifyEqual(G, G_expected, 'AbsTol', 1e-9, 'RelTol', 1e-6)

% Check derivatives wrt parameters using finite differences
sim = SimulateSystem(m, con, 0, opts);
int = sim.int;
int.UseParams = opts.UseParams; % Add UseParams field
t = 0;
verifyDerivativesParameters(a, obj, int)
end

function testObjectiveLogNormalPriorOnSeedParameters(a)
[m, con, ~, opts] = simple_model(); % objectiveWeightedSumOfSquares is the default obj fun

ns = m.ns;
nTs = sum(opts.UseSeeds);
sbar = rand(ns,1);
Vlogbar = rand(ns);
Vlogbar = Vlogbar*Vlogbar.'; % Ensure symmetric and positive definite
Vbar = diag(sbar)*Vlogbar*diag(sbar); % Multiply by square parameters because objective function takes "non-normalized" inputs

obj = objectiveLogNormalPriorOnSeedParameters(sbar, Vbar, 'Testing log seed priors');

% Check objective value
G = ObjectiveValue(m, con, obj, opts);
deltaTs = log(m.s)-log(sbar); 
deltaTs = deltaTs(opts.UseSeeds);
G_expected = deltaTs.'/Vlogbar(opts.UseSeeds,opts.UseSeeds)*deltaTs;
a.verifyEqual(G, G_expected, 'AbsTol', 1e-9, 'RelTol', 1e-6)

% Check derivatives wrt parameters using finite differences
sim = SimulateSystem(m, con, 0, opts);
int = sim.int;
int.UseSeeds = opts.UseSeeds; % Add UseParams field
t = 0;
verifyDerivativesSeeds(a, obj, int)
end

function verifyDerivatives(a, obj, int, t)
ind = find(int.t == t);
x0 = int.y(:,ind);

S.type = '()';
S.subs = {':',ind};
f = @(y)obj.G(setfield(int, 'y', subsasgn(int.y, S, y)));
dfdx = @(y)obj.dGdy(t, setfield(int, 'y', subsasgn(int.y, S, y)));

verifyClose(a, x0, f, dfdx)

S.type = '()';
S.subs = {':',ind};
f = @(y)obj.dGdy(t, setfield(int, 'y', subsasgn(int.y, S, y)));
dfdx = @(y)obj.d2Gdy2(t, setfield(int, 'y', subsasgn(int.y, S, y)));

verifyClose(a, x0, f, dfdx)
end

function verifyDerivativesParameters(a, obj, int)
x0 = int.k;

S.type = '()';
S.subs = {':',1};
f = @(k)obj.G(setfield(int, 'k', subsasgn(int.k, S, k)));
dfdx = @(k)obj.dGdk(setfield(int, 'k', subsasgn(int.k, S, k)));

verifyClose(a, x0, f, dfdx)

S.type = '()';
S.subs = {':',1};
f = @(k)obj.dGdk(setfield(int, 'k', subsasgn(int.k, S, k)));
dfdx = @(k)obj.d2Gdk2(setfield(int, 'k', subsasgn(int.k, S, k)));

verifyClose(a, x0, f, dfdx)
end

function verifyDerivativesSeeds(a, obj, int)
x0 = int.s;

S.type = '()';
S.subs = {':',1};
f = @(s)obj.G(setfield(int, 's', subsasgn(int.s, S, s)));
dfdx = @(s)obj.dGds(setfield(int, 's', subsasgn(int.s, S, s)));

verifyClose(a, x0, f, dfdx)

S.type = '()';
S.subs = {':',1};
f = @(s)obj.dGds(setfield(int, 's', subsasgn(int.s, S, s)));
dfdx = @(s)obj.d2Gds2(setfield(int, 's', subsasgn(int.s, S, s)));

verifyClose(a, x0, f, dfdx)
end

function verifyClose(a, x0, f, dfdx)
% dfdx_finite returned as a sparse matrix
[~, dfdx_finite, dfdx_analytic] = fdiff(x0, f, dfdx);
a.verifyEqual(sparse(dfdx_finite), sparse(dfdx_analytic), 'RelTol', 0.001, 'AbsTol', 1e-4)
end
