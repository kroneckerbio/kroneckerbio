function sol = integrateSens(m, con, opts)

% Constants
nx = m.nx;
nTk = sum(opts.UseParams);
nTs = sum(opts.UseSeeds);
nTq = sum(opts.UseControls);
nT  = nTk + nTs + nTq;

% Construct system
[der, jac] = constructSystem();

if opts.UseModelSeeds
    s = m.s;
else
    s = con.s;
end

if ~con.SteadyState
    % Initial conditions [x0; vec(dxdT0)]
    x0 = m.dx0ds * s + m.x0c;
    
    % Initial effect of rates on sensitivities is 0
    dxdTk0 = zeros(nx, nTk); % Active rate parameters
    
    % Initial effect of seeds on states is dx0ds
    dxdTx0 = m.dx0ds(:,opts.UseSeeds);
    
    % Initial effect of qs on sensitivities is 0
    dxdTq = zeros(nx, nTq);
    
    % Combine them into a vector
    ic = [x0; vec([dxdTk0, dxdTx0, dxdTq])];
else
    % Run to steady-state first
    ic = steadystateSens(m, con, opts);
end

% Input
if opts.UseModelInputs
    u = m.u;
    q = m.q;
else
    u = con.u;
    q = con.q;
end

% Integrate [x; G; dxdT; dGdv] with respect to time
sol = accumulateOde(der, jac, 0, con.tF, ic, u, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+nx*nT));
sol.u = u;
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;
sol.k = m.k;
sol.s = s;
sol.q = q;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x and dxdT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac] = constructSystem()
        
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
        
        % Derivative of [x; dxdT] with respect to time
        function val = derivative(t, joint, u)
            u = u(t);            
            x = joint(1:nx);
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % xT_ --> x_T
            
            % Compute derivative of dxdT
            dxdTdot = dfdx(t,x,u) * dxdT + dfdT(t,x,u); % f_x * x_T + f_T --> f_T
            
            % Combine
            val = [f(t,x,u); vec(dxdTdot)];
        end
        
        % Jacobian of [x; dxdT] derivative
        function val = jacobian(t, joint, u)
            u = u(t);            
            x = joint(1:nx); % x_
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % x_T
            
            % Compute d/dx(dfdT)
            d2xdxdT = sparse(d2fdx2(t,x,u) * dxdT) + d2fdTdx(t,x,u); % fx_T
            d2xdxdT = spermute132(d2xdxdT, [nx,nx,nT], [nx*nT,nx]);
            
            % Combine
            val = [dfdx(t,x,u), sparse(nx,nx*nT);
                       d2xdxdT, kron(IT, dfdx(t,x,u))];
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
    end
end
