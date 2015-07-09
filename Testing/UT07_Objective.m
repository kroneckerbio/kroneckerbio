function tests = UT07_Objective()
tests = functiontests(localfunctions);
end

function testObjectiveValueSimple(a)
[m, con, obj, opts] = simple_model();

G = ObjectiveValue(m, con, obj, opts);
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

function verifyDerivatives(a, obj, int, t)
x0 = int.y(:,int.t == t);

S.type = '()';
S.subs = {':',1};
f = @(y)obj.G(setfield(int, 'y', subsasgn(int.y, S, y)));
dfdx = @(y)obj.dGdy(t, setfield(int, 'y', subsasgn(int.y, S, y)));

verifyClose(a, x0, f, dfdx)

S.type = '()';
S.subs = {':',1};
f = @(y)obj.dGdy(t, setfield(int, 'y', subsasgn(int.y, S, y)));
dfdx = @(y)obj.d2Gdy2(t, setfield(int, 'y', subsasgn(int.y, S, y)));

verifyClose(a, x0, f, dfdx)
end

function verifyClose(a, x0, f, dfdx)
% dfdx_finite returned as a sparse matrix
[~, dfdx_finite, dfdx_analytic] = fdiff(x0, f, dfdx);
a.verifyEqual(sparse(dfdx_finite), sparse(dfdx_analytic), 'RelTol', 0.001)
end
