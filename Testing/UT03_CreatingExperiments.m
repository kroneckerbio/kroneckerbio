function tests = UT03_CreatingExperiments()
tests = functiontests(localfunctions);
end

function testInitialValueExperimentBlank(a)
m = simple_model();
con = InitialValueExperiment(m);
a.verifyEqual(con.s, m.s);
a.verifyEqual(con.u(1), m.u);
a.verifyEqual(con.d(1), zeros(m.ns,1));
a.verifyEqual(con.Periodic, false);
a.verifyEqual(con.SteadyState, false);
end

function testInitialValueExperimentSeed(a)
m = simple_model();
s = rand(m.ns,1);
con = InitialValueExperiment(m, s);
a.verifyEqual(con.s, s);
a.verifyEqual(con.u(1), m.u);
a.verifyEqual(con.d(1), zeros(m.ns,1));
end

function testExperimentInputConstant(a)
m = simple_model();
u = rand(m.nu,1);
con = InitialValueExperiment(m, [], u);
a.verifyEqual(con.s, m.s);
a.verifyEqual(con.u(5), u);
a.verifyEqual(con.d(1), zeros(m.ns,1));
end

function testExperimentDoseConstant(a)
m = simple_model();
d = rand(m.ns,1);
dos = doseConstant(m, d, 1:10);
con = InitialValueExperiment(m, [], [], dos);
a.verifyEqual(con.s, m.s);
a.verifyEqual(con.u(1), m.u);
a.verifyEqual(con.d(1), d)
a.verifyEqual(con.d(1.5), zeros(m.ns,1))
end

function testSimpleExperiment(a)
[unused, con] = simple_model();
verifyDerivatives(a, con)
end

function verifyDerivatives(a, con)
q0 = rand(con.nq,1) + 1;
h0 = rand(con.nh,1) + 1;
con = con.Update(con.s, q0, h0);

verifyClose(a, q0, @(q)q_wrapper(q, 'u', 0), @(q)q_wrapper(q, 'dudq', 0))
verifyClose(a, h0, @(h)h_wrapper(h, 'd', 0), @(h)h_wrapper(h, 'dddh', 0))
verifyClose(a, q0, @(q)q_wrapper(q, 'dudq', 0), @(q)q_wrapper(q, 'd2udq2', 0))
verifyClose(a, h0, @(h)h_wrapper(h, 'dddh', 0), @(h)h_wrapper(h, 'd2ddh2', 0))

verifyClose(a, q0, @(q)q_wrapper(q, 'u', 4), @(q)q_wrapper(q, 'dudq', 4))
verifyClose(a, h0, @(h)h_wrapper(h, 'd', 4), @(h)h_wrapper(h, 'dddh', 4))
verifyClose(a, q0, @(q)q_wrapper(q, 'dudq', 4), @(q)q_wrapper(q, 'd2udq2', 4))
verifyClose(a, h0, @(h)h_wrapper(h, 'dddh', 4), @(h)h_wrapper(h, 'd2ddh2', 4))


    function val = q_wrapper(q, member, t)
        con_temp = con.Update(con.s, q, h0);
        func = con_temp.(member);
        val = func(t);
    end

    function val = h_wrapper(h, member, t)
        con_temp = con.Update(con.s, q0, h);
        func = con_temp.(member);
        val = func(t);
    end
end

function verifyClose(a, x0, f, dfdx)
[unused, dfdx_finite, dfdx_analytic] = fdiff(x0, f, dfdx);
a.verifyEqual(sparse(dfdx_finite), dfdx_analytic, 'RelTol', 0.001)
end
