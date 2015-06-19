function int = integrateSensComp(m, con, tF, eve, fin, opts)
% Constants
nx = m.nx;
nu = m.nu;
ny = m.ny;
nTk = sum(opts.UseParams);
nTs = sum(opts.UseSeeds);
nTq = sum(opts.UseInputControls);
nTh = sum(opts.UseDoseControls);
nT  = nTk + nTs + nTq + nTh;

dxdTStart = nx+1;
dxdTEnd   = nx+nx*nT;

y = m.y;
u = con.u;
dudq = con.dudq;
dydx = m.dydx;
dydu = m.dydu;
dydk = m.dydk;

dkdT = sparse(find(opts.UseParams),1:nTk,1,m.nk,nT);

% Construct system
[der, jac, del] = constructSystem();

if ~con.SteadyState
    % Initial conditions
    order = 1;
    ic = extractICs(m,con,opts,order);
else
    % Run to steady-state first
    ic = steadystateSens(m, con, opts);
end

% Integrate [f; dfdT] over time
sol = accumulateOdeFwdComp(der, jac, 0, tF, ic, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+nx*nT), del, eve, fin);

% Work down
int.Type = 'Integration.Sensitivity.Complex';
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
int.dydk = m.dydk;

int.nT = nT;
int.UseParams = opts.UseParams;
int.UseSeeds = opts.UseSeeds;
int.UseInputControls = opts.UseInputControls;
int.UseDoseControls = opts.UseDoseControls;

int.t = sol.x;
int.x = @(t)evaluate_state(sol, t);
int.u = con.u;
int.y = @(t)evaluate_output(sol, t);

int.dxdT = @(t)evaluate_state_sensitivity(sol, t);
int.dudT = @evaluate_input_sensitivity;
int.dydT = @(t)evaluate_output_sensitivity(sol, t);

nte = numel(sol.ie);
int.ie = sol.ie;
int.te = sol.xe;
int.xe = sol.ye(1:nx,:);
int.ue = u(int.te);
int.ye = y(int.te, int.xe, int.ue);

int.dxedT = sol.ye(dxdTStart:dxdTEnd,:);
int.duedT = zeros(nu*nT,nte);
int.dyedT = zeros(ny*nT,nte);
for it = 1:nte
    duedq_i = dudq(int.te(it)); % u_q
    duedq_i = duedq_i(:,opts.UseInputControls); % u_Q
    duedT_i = [sparse(nu,nTk+nTs), reshape(duedq_i, nu,nTq), sparse(nu,nTh)]; % u_Q -> u_T
    int.duedT(:,it) = vec(duedT_i); % u_T -> uT_

    dyedx_i = dydx(int.te(it), int.xe(:,it), int.ue(:,it)); % y_x
    dyedu_i = dydu(int.te(it), int.xe(:,it), int.ue(:,it)); % y_u
    dyedk_i = dydk(int.te(it), int.xe(:,it), int.ue(:,it)); % y_k
    dxedT_i = reshape(int.dxedT(:,it), nx,nT); % xT_ -> x_T
    int.dyedT(:,it) = vec(dyedx_i * dxedT_i + dyedu_i * duedT_i + dyedk_i * dkdT); % y_x * x_T + y_u * u_T * y_k * k_T -> y_T -> yT_
end

int.sol = sol;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f and dfdT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        
        IT = speye(nT);

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
        dx0dd   = m.dx0ds;
        x0      = m.x0;
        nd      = m.ns;
        
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
            d_i = d(t);
            dx0dd_i = dx0dd(d_i);
            
            deltax = x0(d_i) - x0(zeros(nd,1));
            
            ddeltaxdh = dx0dd_i * dddh(t);
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

    function val = evaluate_state(sol, t)
        val = devals(sol, t, 1:nx);
    end

    function val = evaluate_output(sol, t)
        val = y(t, devals(sol, t, 1:nx), u(t));
    end

    function val = evaluate_state_sensitivity(sol, t)
        val = devals(sol, t, dxdTStart:dxdTEnd);
    end

    function val = evaluate_input_sensitivity(t)
        nt = numel(t);
        
        val = zeros(nu*nT,nt);
        for i = 1:nt
            dudq_i = dudq(t(i)); % u_q
            dudq_i = dudq_i(:,opts.UseInputControls); % u_Q
            val(:,i) = vec([sparse(nu,nTk+nTs), reshape(dudq_i, nu,nTq), sparse(nu,nTh)]); % u_Q -> u_T -> uT_
        end
    end

    function val = evaluate_output_sensitivity(sol, t)
        nt = numel(t);
        
        val = zeros(ny*nT,nt);
        for i = 1:nt
            x_i = deval(sol, t(i), 1:nx); % x_
            u_i =  u(t(i)); % u_
            dxdT_i = reshape(deval(sol, t(i), dxdTStart:dxdTEnd), nx,nT); % xT_ -> x_T
            dudT_i = reshape(evaluate_input_sensitivity(t(i)), nu,nT); % uT_ -> u_T
            val(:,i) = vec(dydx(t(i), x_i, u_i) * dxdT_i + dydu(t(i), x_i, u_i) * dudT_i + dydk(t(i), x_i, u_i) * dkdT); % y_x * x_T + y_u * u_T * y_k * k_T -> y_T -> yT_
        end
    end
end
