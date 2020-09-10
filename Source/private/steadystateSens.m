function ic = steadystateSens(m, con, opts)

% Constants
nx = m.nx;
nTk = sum(opts.UseParams);
nTs = sum(opts.UseSeeds);
nTq = sum(opts.UseInputControls);
nTh = sum(opts.UseDoseControls);
nT  = nTk + nTs + nTq + nTh;

% Construct system
[der, jac] = constructSystem();

% Initial conditions [x0; vec(dxdT0)]
order = 1;
ic = extractICs(m,con,opts,order);

abstol = opts.AbsTol(1:nx+nx*nT);
reltol = opts.RelTol;
basal_discontinuities = con.private.BasalDiscontinuities;
timescale = con.private.TimeScale;
ic = iterate_steady_state(der, jac, ic, nx, abstol, reltol, basal_discontinuities, timescale);

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f and dfdT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, eve] = constructSystem()
        
        IT = speye(nT);

        dxdTStart = nx+1;
        dxdTEnd   = nx+nx*nT;
        
        f       = m.f;
        dfdx    = m.dfdx;
        dfdu    = m.dfdu;
        dfdk    = m.dfdk;
        dfdT    = @dfdTSub;
        d2fdx2  = m.d2fdx2;
        d2fdudx = m.d2fdudx;
        d2fdkdx = m.d2fdkdx;
        d2fdTdx = @d2fdTdxSub;
        uf      = con.private.basal_u;
        dudq    = con.private.basal_dudq;
        
        der = @derivative;
        jac = @jacobian;
        eve = @events;
        
        % Derivative of [x; dxdT] with respect to time
        function val = derivative(t, joint)
            u = uf(t);            
            x = joint(1:nx);
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % xT_ --> x_T
            
            % Derivative of dxdT
            dxdTdot = dfdx(-1,x,u) * dxdT + dfdT(t,x,u); % f_x * x_T + f_T --> f_T
            
            % Combine
            val = [f(-1,x,u); vec(dxdTdot)];
        end
        
        % Jacobian of [x; dxdT] derivative
        function val = jacobian(t, joint)
            u = uf(t);            
            x = joint(1:nx); % x_
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % x_T
            
            % Compute d/dx(dfdT)
            d2xdxdT = sparse(d2fdx2(-1,x,u) * dxdT) + d2fdTdx(t,x,u); % fx_T
            d2xdxdT = spermute132(d2xdxdT, [nx,nx,nT], [nx*nT,nx]);
            
            % Combine
            val = [dfdx(-1,x,u), sparse(nx,nx*nT);
                        d2xdxdT, kron(IT, dfdx(-1,x,u))];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = dfdk(-1,x,u);
            dfdq = dfdu(-1,x,u) * dudq(t);
            val = [val(:,opts.UseParams), sparse(nx,nTs), dfdq(:,opts.UseInputControls), sparse(nx,nTh)];
        end
        
        % Modifies d2fdkdx to relate only to the parameters of interest
        function val = d2fdTdxSub(t, x, u)
            val = d2fdkdx(-1,x,u); % fx_k
            d2fdqdx = d2fdudx(-1,x,u) * dudq(t); % fx_u * u_q -> fx_q
            val = [val(:,opts.UseParams), sparse(nx*nx, nTs), d2fdqdx(:,opts.UseInputControls), sparse(nx*nx,nTh)]; % fx_T
        end
        
    end
end
