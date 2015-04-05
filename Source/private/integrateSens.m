function sol = integrateSens(m, con, opts)

% Constants
nx = m.nx;
nTk = sum(opts.UseParams);
nTs = sum(opts.UseSeeds);
nTq = sum(opts.UseInputControls);
nTh = sum(opts.UseDoseControls);
nT  = nTk + nTs + nTq + nTh;

% Construct system
[der, jac, del] = constructSystem();

if ~con.SteadyState
    % Initial conditions [x0; vec(dxdT0)]
    x0 = m.dx0ds * con.s + m.x0c;
    
    % Initial effect of rates on sensitivities is 0
    dxdTk = zeros(nx, nTk); % Active rate parameters
    
    % Initial effect of seeds on states is dx0ds
    dxdTs = m.dx0ds(:,opts.UseSeeds);
    
    % Initial effect of qs on sensitivities is 0
    dxdTq = zeros(nx, nTq);
    
    % Initial effect of hs on sensitivities is 0
    dxdTh = zeros(nx, nTh);

    % Combine them into a vector
    ic = [x0; vec([dxdTk, dxdTs, dxdTq, dxdTh])];
else
    % Run to steady-state first
    ic = steadystateSens(m, con, opts);
end

% Integrate [f; dfdT] over time
sol = accumulateOdeFwd(der, jac, 0, con.tF, ic, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+nx*nT), del);
sol.nx = nx;
sol.nT = nT;
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f and dfdT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        
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
        uf      = con.u;
        dudq    = con.dudq;
        d       = con.d;
        dddh    = con.dddh;
        dx0ds   = m.dx0ds;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [x; dxdT] with respect to time
        function val = derivative(t, joint)
            u = uf(t);            
            x = joint(1:nx);
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % xT_ -> x_T
            
            % Derivative of dxdT
            dxdTdot = dfdx(t,x,u) * dxdT + dfdT(t,x,u); % f_x * x_T + f_T -> f_T
            
            % Combine
            val = [f(t,x,u); vec(dxdTdot)];
        end
        
        % Jacobian of [x; dxdT] derivative
        function val = jacobian(t, joint)
            u = uf(t);            
            x = joint(1:nx); % x_
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % x_T
            
            % Compute d/dx(dfdT)
            d2fdxdT = sparse(d2fdx2(t,x,u) * dxdT) + d2fdTdx(t,x,u); % fx_T
            d2fdxdT = spermute132(d2fdxdT, [nx,nx,nT], [nx*nT,nx]);
            
            % Combine
            val = [dfdx(t,x,u), sparse(nx,nx*nT);
                       d2fdxdT, kron(IT, dfdx(t,x,u))];
        end
        
        % Dosing
        function val = delta(t, joint)
            deltax = dx0ds * d(t);
            
            ddeltaxdh = dx0ds * dddh(t);
            ddeltaxdT = [zeros(nx,nTk+nTs+nTq), ddeltaxdh(:,opts.UseDoseControls)];
            
            val = [deltax; vec(ddeltaxdT)];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = m.dfdk(t,x,u);
            dfdq = dfdu(t,x,u) * dudq(t);
            val = [val(:,opts.UseParams), sparse(nx,nTs), dfdq(:,opts.UseInputControls), sparse(nx,nTh)];
        end
        
        % Modifies d2fdkdx to relate only to the parameters of interest
        function val = d2fdTdxSub(t, x, u)
            val = m.d2fdkdx(t,x,u); % fx_k
            d2fdqdx = d2fdudx(t,x,u) * dudq(t); % fx_u * u_q -> fx_q
            val = [val(:,opts.UseParams), sparse(nx*nx, nTs), d2fdqdx(:,opts.UseInputControls), sparse(nx*nx,nTh)]; % fx_T
        end
    end
end
