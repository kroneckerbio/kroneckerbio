function ic = steadystateSys(m, con, opts)

% Constants
nx = m.nx;

% Construct system
[der, jac] = constructSystem();

order = 0;
ic = extractICs(m,con,opts,order);

abstol = opts.AbsTol(1:nx);
reltol = opts.RelTol;
basal_discontinuities = con.private.BasalDiscontinuities;
timescale = con.private.TimeScale;
ic = iterate_steady_state(der, jac, ic, nx, abstol, reltol, basal_discontinuities, timescale);

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, eve] = constructSystem()
        
        f    = m.f;
        dfdx = m.dfdx;
        uf    = con.private.basal_u;
        
        der = @derivative;
        jac = @jacobian;
        eve = @events;
        
        % Derivative of x with respect to time
        function val = derivative(t, x)
            u   = uf(t);
            val = f(-1, x, u);
        end
        
        % Jacobian of x derivative
        function val = jacobian(t, x)
            u   = uf(t);
            val = dfdx(-1, x, u);
        end
        
    end

end