function sol = integrateSysSelect(m, con, tGet, opts)
% Constants
nx = m.nx;

% Construct system
[der, jac] = constructSystem();

% Initial conditions
if ~con.SteadyState
    if opts.UseModelICs
        ic = m.x0;
    else
        ic = con.x0;
    end
else
    ic = steadystateSys(m, con, opts);
end

% Input
if opts.UseModelInputs
    u = m.u;
else
    u = con.u;
end

% Integrate x over time
sol = accumulateOde(der, jac, 0, con.tF, ic, u, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx), [], 1, [], [], [], tGet);
sol.u = u(tGet);
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac] = constructSystem()
        f    = m.f;
        dfdx = m.dfdx;
        
        der = @derivative;
        jac = @jacobian;

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
    end
end