function sol = integrateSysSelect(m, con, tGet, opts)

% Constants
nx = m.nx;

% Construct system
[der, jac, del] = constructSystem();

% Initial conditions
if opts.UseModelSeeds
    s = m.s;
else
    s = con.s;
end

if ~con.SteadyState
    ic = m.dx0ds * s + m.x0c;
else
    ic = steadystateSys(m, con, opts);
end

% Input
if opts.UseModelInputs
    u = m.u;
    q = m.q;
else
    u = con.u;
    q = con.q;
end

% Integrate x over time
sol = accumulateOdeFwdSelect(der, jac, 0, con.tF, ic, u, con.Discontinuities, tGet, 1:nx, opts.RelTol, opts.AbsTol(1:nx), del);
sol.u = u(tGet);
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;
sol.k = m.k;
sol.s = s;
sol.q = q;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        f    = m.f;
        dfdx = m.dfdx;
        d    = con.d;
        dx0ds = m.dx0ds;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;

        % Derivative of x with respect to time
        function val = derivative(t, x, u)
            u   = u(t);
            val = f(t, x, u);
        end

        % Jacobian of x derivative
        function val = jacobian(t, x, u)
            u   = u(t);
            val = dfdx(t, x, u);
        end
        
        % Dosing
        function val = delta(t, x)
            val = dx0ds * d(t);
        end
    end
end