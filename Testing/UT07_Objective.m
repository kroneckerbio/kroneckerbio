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
verifyDerivatives(a, m, obj, int, t)
end

function testObjectiveWeightedSumOfSquaresNonNeg(a)
[m, con, obj, opts] = simple_model('objectiveWeightedSumOfSquaresNonNeg');

obs = observationSelect(1:6);
sim = SimulateSystem(m, con, obs, opts);

int = sim.int;
t = 1;
verifyDerivatives(a, m, obj, int, t)
end

function verifyDerivatives(a, m, obj, int, t)
x0 = int.y(:,int.t == t);

f = @(x)G_wrapper(int, x);
dfdx = @(x)dGdx_wrapper(int, x);
verifyClose(a, x0, f, dfdx)

    function val = G_wrapper(int_i, x)
        u = int_i.u(:,int_i.t == t);
        int_i.x(:,int_i.t == t) = x;
        int_i.y(:,int_i.t == t) = m.y(t, x, u);
        val = obj.G(int_i);
    end

    function val = dGdx_wrapper(int_i, x)
        u = int_i.u(:,int_i.t == t);
        int_i.x(:,int_i.t == t) = x;
        int_i.y(:,int_i.t == t) = m.y(t, x, u);
        val = obj.dGdx(t, int_i);
    end
end

function verifyClose(a, x0, f, dfdx)
% dfdx_finite returned as a sparse matrix
[~, dfdx_finite, dfdx_analytic] = fdiff(x0, f, dfdx);
a.verifyEqual(full(dfdx_finite), dfdx_analytic, 'RelTol', 0.001)
end
