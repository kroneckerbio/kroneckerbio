function sol = integrateSys(m, con, opts)

% Constants
nx = m.nx;

% Construct system
[der, jac, del] = constructSystem();

if ~con.SteadyState
    ic = m.dx0ds * con.s + m.x0c;
else
    ic = steadystateSys(m, con, opts);
end

% Integrate f over time
sol = accumulateOdeFwd(der, jac, 0, con.tF, ic, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx), del);
sol.nx = nx;
sol.u = con.u;
sol.y_ = m.y;
sol.dydx = m.dydx;
sol.dydu = m.dydu;
sol.k = m.k;
sol.s = con.s;
sol.q = con.q;
sol.h = con.h;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        f     = m.f;
        dfdx  = m.dfdx;
        uf    = con.u;
        d     = con.d;
        dx0ds = m.dx0ds;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of x with respect to time
        function val = derivative(t, x)
            u   = uf(t);
            val = f(t, x, u);
        end
        
        % Jacobian of x derivative
        function val = jacobian(t, x)
            u   = uf(t);
            val = dfdx(t, x, u);
        end
        
        % Dosing
        function val = delta(t, x)
            val = dx0ds * d(t);
        end
    end
end
