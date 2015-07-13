function ic = steadystateSys(m, con, opts)

% Constants
nx = m.nx;

% Construct system
[der, jac, eve] = constructSystem();

order = 0;
ic = extractICs(m,con,opts,order);

% Integrate f over time
%     accumulateOdeFwdSimp(der, jac, t0, tF, ic, discontinuities, t_get, nonnegative, RelTol, AbsTol, delta, events, is_finished)
sol = accumulateOdeFwdSimp(der, jac, 0, inf, ic, con.private.BasalDiscontinuities, 0, 1:nx, opts.RelTol, opts.AbsTol(1:nx), [], eve, @(cum_sol)true);

% Return steady-state value
ic = sol.ye;

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
        
        % Steady-state event
        function [value, isTerminal, direction] = events(t, x)
            u = uf(t);

            % Absolute change
            absDiff = con.private.TimeScale * f(-1,x,u); % Change over an entire simulation
            
            % Relative change
            relDiff = absDiff ./ x;
            
            % Either absolute change or relative change must be less than
            % the tolerance for all species
            value = max(min(abs(absDiff) - opts.AbsTol(1:nx), abs(relDiff) - opts.RelTol));
            if value < 0
                value = 0;
            end
            
            % Always end and only care about drops below the threshold
            isTerminal = true;
            direction = -1;
        end
    end
end
