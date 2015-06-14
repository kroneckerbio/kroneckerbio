function int = integrateSysComp(m, con, tF, eve, fin, opts)

% Constants
nx = m.nx;

y = m.y;
u = con.u;

% Construct system
[der, jac, del] = constructSystem();

if ~con.SteadyState
    order = 0;
    ic = extractICs(m,con,opts,order);
else
    ic = steadystateSys(m, con, opts);
end

% Integrate f over time
sol = accumulateOdeFwdComp(der, jac, 0, tF, ic, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx), del, eve, fin);

% Work down
int.Type = 'Integration.System.Complex';
int.Name = [m.Name ' in ' con.Name];

int.x_names = vec({m.States.Name});
int.u_names = vec({m.Inputs.Name});
int.y_names = vec({m.Outputs.Name});
int.k_names = vec({m.Parameters.Name});
int.s_names = vec({m.Seeds.Name});

int.nx = nx;
int.ny = m.ny;
int.nu = m.nu;
int.nk = m.nk;
int.ns = m.ns;
int.nq = con.nq;
int.nh = con.nh;
int.k = m.k;
int.s = con.s;
int.q = con.q;
int.h = con.h;

int.dydx = m.dydx;
int.dydu = m.dydu;

int.t = sol.x;
int.x = @(t)devals(sol, t);
int.u = con.u;
int.y = @(t)y(t, devals(sol, t), u(t));

int.ie = sol.ie;
int.te = sol.xe;
int.xe = sol.ye;
int.ue = u(int.te);
int.ye = y(int.te, int.xe, int.ue);

int.sol = sol;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        f     = m.f;
        dfdx  = m.dfdx;
        u     = con.u;
        d     = con.d;
        x0    = m.x0;
        nd    = m.ns;
        
        y = m.y;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of x with respect to time
        function val = derivative(t, x)
            u_t = u(t);
            val = f(t, x, u_t);
        end
        
        % Jacobian of x derivative
        function val = jacobian(t, x)
            u_t = u(t);
            val = dfdx(t, x, u_t);
        end
        
        % Dosing
        function val = delta(t, x)
            val = x0(d(t)) - x0(zeros(nd,1));
        end
    end
end
