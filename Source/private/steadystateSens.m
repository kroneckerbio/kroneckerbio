function ic = steadystateSens(m, con, opts)

% Constants
nx = m.nx;
nTk = sum(opts.UseParams);
nTs = sum(opts.UseSeeds);
nTq = sum(opts.UseControls);
nT  = nTk + nTs + nTq;

% Construct system
[der, jac, eve] = constructSystem();

% Initial conditions [x0; vec(dxdT0)]
if opts.UseModelSeeds
    x0 = m.dx0ds * m.s + m.x0c;
else
    x0 = m.dx0ds * con.s + m.x0c;
end

% Initial effect of rates on sensitivities is 0
dxdTk = zeros(nx, nTk); % Active rate parameters

% Initial effect of seeds on states is dx0ds
dxdTs = m.dx0ds(:,opts.UseSeeds);

% Initial effect of qs on sensitivities is 0
dxdTq = zeros(nx, nTq);

% Combine them into a vector
ic = [x0; vec([dxdTk, dxdTs, dxdTq])];

% Input
if opts.UseModelInputs
    u = m.u;
else
    u = con.u;
end

% Integrate [x; dxdT] with respect to time
sol = accumulateOde(der, jac, 0, inf, ic, u, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+nx*nT), [], 1, eve, [], 1, 0);

% Return steady-state value
ic = sol.ye;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x and dxdT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, eve] = constructSystem()
        
        IT = speye(nT);

        dxdTStart = nx+1;
        dxdTEnd   = nx+nx*nT;
        
        f       = m.f;
        dfdx    = m.dfdx;
        dfdu    = m.dfdu;
        dfdT    = @dfdTSub;
        d2fdx2  = m.d2fdx2;
        d2fdudx = m.d2fdudx;
        d2fdTdx = @d2fdTdxSub;
        if opts.UseModelInputs
            dudq = m.dudq;
        else
            dudq = con.dudq;
        end
        
        der = @derivative;
        jac = @jacobian;
        eve = @events;
        
        % Derivative of [x; dxdT] with respect to time
        function val = derivative(t, joint, u)
            u = u(-1);            
            x = joint(1:nx);
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % xT_ --> x_T
            
            % Compute derivative of dxdT
            dxdTdot = dfdx(-1,x,u) * dxdT + dfdT(-1,x,u); % f_x * x_T + f_T --> f_T
            
            % Combine
            val = [f(-1,x,u); vec(dxdTdot)];
        end
        
        % Jacobian of [x; dxdT] derivative
        function val = jacobian(t, joint, u)
            u = u(-1);            
            x = joint(1:nx); % x_
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % x_T
            
            % Compute d/dx(dfdT)
            d2xdxdT = sparse(d2fdx2(-1,x,u) * dxdT) + d2fdTdx(-1,x,u); % fx_T
            d2xdxdT = spermute132(d2xdxdT, [nx,nx,nT], [nx*nT,nx]);
            
            % Combine
            val = [dfdx(-1,x,u), sparse(nx,nx*nT);
                        d2xdxdT, kron(IT, dfdx(-1,x,u))];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = m.dfdk(t,x,u);
            dfdq = dfdu(t,x,u) * dudq(t);
            val = [val(:,opts.UseParams), sparse(nx, nTs), dfdq(:,opts.UseControls)];
        end
        
        % Modifies d2fdkdx to relate only to the parameters of interest
        function val = d2fdTdxSub(t, x, u)
            val = m.d2fdkdx(t,x,u);
            d2fdqdx = d2fdudx(t,x,u) * dudq(t);
            val = [val(:,opts.UseParams), sparse(nx*nx, nTs), d2fdqdx(:,opts.UseControls)];
        end
        
        % Steady-state event
        function [value, isTerminal, direction] = events(t, joint, u)
            u = u(-1);
            x = joint(1:nx); % x_

            % Absolute change
            absDiff = con.tF * f(-1,x,u); % Change over an entire simulation
            
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