function tests = UT07_Objective()
tests = functiontests(localfunctions);
end

function testObjectiveValueSimple(a)
[m, con, obj, opts] = simple_model();

G = ObjectiveValue(m, con, obj, opts);
end

function testObjectiveWeightedSumOfSquares(a)
[m, con, obj, opts] = simple_model();

sim = SimulateSelect(m, con, 1:6, opts);

verifyDerivatives(a, obj, sim.sol, 1)

end

function verifyDerivatives(a, obj, sol, t)
x0 = sol.y(:,sol.x == t);
rat=1; fin=1; ana=1;
verifyClose(a, x0, @(x)G_wrapper(sol, x), @(x)dGdx_wrapper(sol, x))

    function val = G_wrapper(sol_i, x)
        sol_i.y(:,sol_i.x == t) = x;
        val = obj.G(sol_i);
    end

    function val = dGdx_wrapper(sol_i, x)
        sol_i.y(:,sol_i.x == t) = x;
        val = obj.dGdx(t, sol_i);
    end
end

function verifyClose(a, x0, f, dfdx)
[unused, dfdx_finite, dfdx_analytic] = fdiff(x0, f, dfdx);
a.verifyEqual(dfdx_finite, dfdx_analytic, 'RelTol', 0.001)
end
