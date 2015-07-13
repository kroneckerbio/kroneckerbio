function tests = UT02_LoadingModels()
tests = functiontests(localfunctions);
end

function testEquilibriumLoading(a)
m = LoadModelMassAction('Equilibrium.txt');
a.verifyEqual(m.Name, 'Equilibrium')
end

function testSimpleLoading(a)
m = LoadModelMassAction('Simple.txt');
a.verifyEqual(m.Name, 'Simple')
verifyDerivatives(a, m)
end

function verifyDerivatives(a, m)
% Make all the values about the same size
x0 = rand(m.nx,1) + 1;
u0 = rand(m.nu,1) + 1;
k0 = rand(m.nk,1) + 1;
m = m.Update(k0);

verifyClose(a, x0, @(x)m.f(0,x,u0), @(x)m.dfdx(0,x,u0))
verifyClose(a, u0, @(u)m.f(0,x0,u), @(u)m.dfdu(0,x0,u))
verifyClose(a, k0, @(k)k_wrapper(k, 'f'), @(k)k_wrapper(k, 'dfdk'))

verifyClose(a, x0, @(x)m.dfdx(0,x,u0), @(x)m.d2fdx2(0,x,u0))
verifyClose(a, u0, @(u)m.dfdu(0,x0,u), @(u)m.d2fdu2(0,x0,u))
verifyClose(a, k0, @(k)k_wrapper(k, 'dfdk'), @(k)k_wrapper(k, 'd2fdk2'))

verifyClose(a, u0, @(u)m.dfdx(0,x0,u), @(u)m.d2fdudx(0,x0,u))
verifyClose(a, x0, @(x)m.dfdu(0,x,u0), @(x)m.d2fdxdu(0,x,u0))
verifyClose(a, k0, @(k)k_wrapper(k, 'dfdx'), @(k)k_wrapper(k, 'd2fdkdx'))

verifyClose(a, x0, @(x)m.dfdk(0,x,u0), @(x)m.d2fdxdk(0,x,u0))
verifyClose(a, k0, @(k)k_wrapper(k, 'dfdu'), @(k)k_wrapper(k, 'd2fdkdu'))
verifyClose(a, u0, @(u)m.dfdk(0,x0,u), @(u)m.d2fdudk(0,x0,u))

    function val = k_wrapper(k, member)
        m_temp = m.Update(k);
        func = m_temp.(member);
        val = func(0, x0, u0);
    end
end

function verifyClose(a, x0, f, dfdx)
[~, dfdx_finite, dfdx_analytic] = fdiff(x0, f, dfdx);
a.verifyEqual(sparse(dfdx_finite), dfdx_analytic, 'RelTol', 0.001)
end
