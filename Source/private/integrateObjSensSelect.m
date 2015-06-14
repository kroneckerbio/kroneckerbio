function sol = integrateObjSensSelect(m, con, obj, tGet, opts)

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
    % Initial conditions
    order = 1;
    x0_dx0dT = extractICs(m,con,opts,order);
    
    % Combine them into a vector
    ic = [x0_dx0dT(1:nx); 0; x0_dx0dT(nx+1:end); zeros(nT,1)];
else
    % Run to steady-state first
    ic = steadystateSens(m, con, opts);
    
    % Include objective and gradient
    ic = [ic(1:nx); 0; ic(nx+1:nx+1+nx*nT); zeros(nT,1)];
end

% Integrate [f; g; dfdT; dgdT] over time
sol = accumulateOde(der, jac, 0, con.tF, ic, con.Discontinuities, tGet, 1:nx, opts.RelTol, opts.AbsTol(1:nx+1+nx*nT+nT), del);
sol.u = con.u(tGet);
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;
sol.k = m.k;
sol.s = con.s;
sol.q = con.q;
sol.h = con.h;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f, g, dfdT, and dgdT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        
        IT = speye(nT);

        dxdTStart = nx+1+1;
        dxdTEnd   = nx+1+nx*nT;
        
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
        dx0dd   = m.dx0ds;
        x0      = m.x0;
        nd      = m.ns;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [x; G; dxdT; dGdT] with respect to time
        function val = derivative(t, joint)
            u = uf(t);            
            x = joint(1:nx);
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % xT_ -> x_T
            
            % Derivative of dxdT
            dxdTdot = dfdx(t,x,u) * dxdT + dfdT(t,x,u); % f_x * x_T + f_T -> f_T
            
            % Sum continuous objective functions
            g = 0;
            dgdT = zeros(nT,1);
            for i = 1:nObj
                g = g + opts.ObjWeights(i) * obj(i).g(t,x,u);
                dgdT = dgdT + opts.ObjWeights(i) * vec(vec(obj(i).dgdx(t,x,u)).' * dxdT); % T_ + (_x * x_T -> _T -> T_) + T_ -> T_
            end
                        
            val = [f(t,x,u); g; vec(dxdTdot); dgdT];
        end
        
        % Jacobian of [x; G; dxdT; dGdv] derivative
        function val = jacobian(t, joint)
            u = uf(t);            
            x = joint(1:nx); % x_
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % x_T
            
            % Sum continuous objective gradients
            dgdx = zeros(1,nx); % _x
            d2gdxdT = sparse(nT,nx); % T_x
            for i = 1:nObj
                dgdx = dgdx + opts.ObjWeights(i) * vec(obj(i).dgdx(t,x,u)).'; % _x
                
                % Compute d/dx(dgdT)
                d2gdxdT = d2gdxdT + opts.ObjWeights(i) * (dxdT.' * obj(i).d2gdx2(t,x,u)); % T_x + T_x * x_x + T_x -> T_x
            end
            
            % Compute d/dx(dfdT)
            d2xdxdT = d2fdx2(t,x,u) * dxdT + d2fdTdx(t,x,u); % fx_T
            d2xdxdT = full(d2xdxdT); % fx_T
            d2xdxdT = reshape(d2xdxdT, nx,nx,nT); % f_x_T
            d2xdxdT = permute(d2xdxdT, [1,3,2]); % f_T_x
            d2xdxdT = reshape(d2xdxdT, nx*nT,nx); % fT_x
            
            % Combine
            val = [dfdx(t,x,u), sparse(nx,1+nx*nT+nT);
                          dgdx, sparse(1,1+nx*nT+nT);
                       d2xdxdT, sparse(nx*nT,1), kron(IT, dfdx(t,x,u)), sparse(nx*nT,nT);
                       d2gdxdT, sparse(nT,1),    kron(IT, dgdx),        sparse(nT,nT)];
        end
        
        % Dosing
        function val = delta(t, joint)
            d_i = d(t);
            dx0dd_i = dx0dd(d_i);
            
            deltax = x0(d_i) - x0(zeros(nd,1));
            
            ddeltaxdh = dx0dd_i * dddh(t);
            ddeltaxdT = [zeros(nx,nTk), zeros(nx,nTs), zeros(nx,nTq), ddeltaxdh(:,opts.UseDoseControls)];
            
            val = [deltax; 0; vec(ddeltaxdT); zeros(nT,1)];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = m.dfdk(t,x,u);
            dfdq = dfdu(t,x,u) * dudq(t);
            val = [val(:,opts.UseParams), sparse(nx,nTs), dfdq(:,opts.UseInputControls), sparse(nx,nTq)];
        end
        
        % Modifies d2fdkdx to relate only to the parameters of interest
        function val = d2fdTdxSub(t, x, u)
            val = m.d2fdkdx(t,x,u);
            d2fdqdx = d2fdudx(t,x,u) * dudq(t);
            val = [val(:,opts.UseParams), sparse(nx*nx, nTs), d2fdqdx(:,opts.UseInputControls), sparse(nx*nx,nTq)];
        end
    end
end
