function int = integrateSensSimp(m, con, tF, eve, fin, t_get, opts)

% Constants
nx = m.nx;
nu = m.nu;
ny = m.ny;
nTk = sum(opts.UseParams);
nTs = sum(opts.UseSeeds);
nTq = sum(opts.UseInputControls);
nTh = sum(opts.UseDoseControls);
nT  = nTk + nTs + nTq + nTh;

y = m.y;
u = con.u;
dudq = con.dudq;
dydx = m.dydx;
dydu = m.dydu;

nt = numel(t_get);

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
sol = accumulateOdeFwdSimp(der, jac, 0, tF, ic, con.Discontinuities, t_get, 1:nx, opts.RelTol, opts.AbsTol(1:nx+nx*nT), del, eve, fin);

% Work down
int.Type = 'Integration.Sensitivity.Simple';
int.Name = [m.Name ' in ' con.Name];

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

int.nT = nT;
int.UseParams = opts.UseParams;
int.UseSeeds = opts.UseSeeds;
int.UseInputControls = opts.UseInputControls;
int.UseDoseControls = opts.UseDoseControls;

int.t = sol.x;
int.x = sol.y(1:nx,:);
int.u = con.u(int.t);
int.y = y(int.t, int.x, int.u);

int.dxdT = sol.y(nx+1:end,:);
int.dudT = zeros(nu*nT,nt);
int.dydT = zeros(ny*nT,nt);
for it = 1:nt
    dudq_i = dudq(int.t(it)); % u_q
    dudq_i = dudq_i(:,opts.UseInputControls); % u_Q
    dudT_i = [zeros(nu,nTk+nTs), reshape(dudq_i, nu,nTq), zeros(nu,nTh)]; % u_Q -> u_T
    int.dudT(:,it) = vec(dudT_i); % u_T -> uT_
    
    dydx_i = dydx(int.t(it), int.x(:,it), int.u(:,it)); % y_x
    dydu_i = dydu(int.t(it), int.x(:,it), int.u(:,it)); % y_u
    dxdT_i = reshape(int.dxdT(:,it), nx,nT); % xT_ -> x_T
    int.dydT(:,it) = vec(dydx_i * dxdT_i + dydu_i * dudT_i); % y_x * x_T + y_u * u_T -> y_T -> yT_
end

nte = numel(sol.ie);
int.ie = sol.ie;
int.te = sol.xe;
int.xe = sol.ye(1:nx,:);
int.ue = u(int.te);
int.ye = y(int.te, int.xe, int.ue);

int.dxedT = sol.ye(nx+1:end,:);
int.duedT = zeros(nu*nT,nte);
int.dyedT = zeros(ny*nT,nte);
for it = 1:nte
    duedq_i = dudq(int.te(it)); % u_q
    duedq_i = duedq_i(:,opts.UseInputControls); % u_Q
    duedT_i = [zeros(nu,nTk+nTs), reshape(duedq_i, nu,nTq), zeros(nu,nTh)]; % u_Q -> u_T
    int.duedT(:,it) = vec(duedT_i); % u_T -> uT_

    dyedx_i = dydx(int.te(it), int.xe(:,it), int.ue(:,it)); % y_x
    dyedu_i = dydu(int.te(it), int.xe(:,it), int.ue(:,it)); % y_u
    dxedT_i = reshape(int.dxedT(:,it), nx,nT); % xT_ -> x_T
    int.dyedT(:,it) = vec(dyedx_i * dxedT_i + dyedu_i * duedT_i); % y_x * x_T + y_u * u_T -> y_T -> yT_
end

int.sol = sol;

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
        dudq    = con.dudq;
        d       = con.d;
        dddh    = con.dddh;
        dx0ds   = m.dx0ds;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [x; dxdT] with respect to time
        function val = derivative(t, joint)
            u_t = u(t);            
            x = joint(1:nx);
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % xT_ -> x_T
            
            % Derivative of dxdT
            dxdTdot = dfdx(t,x,u_t) * dxdT + dfdT(t,x,u_t); % f_x * x_T + f_T -> f_T
            
            % Combine
            val = [f(t,x,u_t); vec(dxdTdot)];
        end
        
        % Jacobian of [x; dxdT] derivative
        function val = jacobian(t, joint)
            u_t = u(t);            
            x = joint(1:nx); % x_
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % x_T
            
            % Compute d/dx(dfdT)
            d2fdxdT = sparse(d2fdx2(t,x,u_t) * dxdT) + d2fdTdx(t,x,u_t); % fx_T
            d2fdxdT = spermute132(d2fdxdT, [nx,nx,nT], [nx*nT,nx]);
            
            % Combine
            val = [dfdx(t,x,u_t), sparse(nx,nx*nT);
                       d2fdxdT, kron(IT, dfdx(t,x,u_t))];
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
