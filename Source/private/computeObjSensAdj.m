function [G, D] = computeObjSensAdj(m, con, obj, opts)
verboseAll = max(opts.Verbose-1,0);
if verboseAll; tic; end

% Constants
nx = m.nx;
nTk = sum(opts.UseParams);
nTs = sum(sum(opts.UseSeeds));
nTq = sum(cat(1,opts.UseInputControls{:}));
nTh = sum(cat(1,opts.UseDoseControls{:}));
nT  = nTk + nTs + nTq +nTh;
n_con = numel(con);
n_obj = size(obj,1);

y = m.y;
dydx = m.dydx;
dydu = m.dydu;
dydk = m.dydk;

% Initialize variables
G = 0;
D = zeros(nT,1);
Tsind = nTk; % Stores the position in D where the first x0 parameter goes for each iCon
Tqind = nTk+nTs; % Stores the position in D where the first q parameter goes for each iCon
Thind = nTk+nTs+nTq; % Stores the position in D where the first h parameter goes for each iCon

if opts.Verbose; disp('Integrating adjoint...'); end
for i_con = 1:n_con
    if verboseAll; tic; end
    opts_i = opts;
    
    % Modify opts structure
    opts_i.AbsTol = opts.AbsTol{i_con};
    opts_i.ObjWeights = opts.ObjWeights(:,i_con);

    UseSeeds_i = opts.UseSeeds(:,i_con);
    opts_i.UseSeeds = UseSeeds_i;
    inTs = nnz(UseSeeds_i);
    
    UseInputControls_i = opts.UseInputControls{i_con};
    opts_i.UseInputControls = UseInputControls_i;
    inTq = nnz(UseInputControls_i);
    
    UseDoseControls_i = opts.UseDoseControls{i_con};
    opts_i.UseDoseControls = UseDoseControls_i;
    inTh = nnz(UseDoseControls_i);
    
    inT = nTk + inTs + inTq + inTh;
    
    % Seeds
    s = con(i_con).s;
    
    % Input
    u = con(i_con).u;
    d = con(i_con).d;
    
    % * Integrate to steady-state
    if con(i_con).SteadyState
        ssSol = integrateSteadystateSys(m, con(i_con), opts_i);
        ic = ssSol.ye(:,end);
    else
        order = 0;
        ic = extractICs(m,con(i_con),opts_i,order);
    end
    
    [tF, eve, fin] = collectObservations(m, con(i_con), obj(:,i_con));

    % * Integrate system *
    % Do not use select methods since the solution is needed at all time
    if opts.continuous(i_con)
        [der, jac, del] = constructObjectiveSystem();
        sol_sys = accumulateOdeFwdComp(der, jac, 0, tF, [ic; 0], con(i_con).Discontinuities, 1:nx, opts.RelTol, opts_i.AbsTol(1:nx+1), del, eve, fin);
    else
        [der, jac, del] = constructSystem();
        sol_sys = accumulateOdeFwdComp(der,jac, 0, tF, ic, con(i_con).Discontinuities, 1:nx, opts.RelTol, opts_i.AbsTol(1:nx), del, eve, fin);
    end
    
    % Work down
    int_sys = struct;
    int_sys.Type = 'Integration.System.Complex';
    int_sys.Name = [m.Name ' in ' con(i_con).Name];

    int_sys.nx = nx;
    int_sys.ny = m.ny;
    int_sys.nu = m.nu;
    int_sys.nk = m.nk;
    int_sys.ns = m.ns;
    int_sys.nq = con(i_con).nq;
    int_sys.nh = con(i_con).nh;
    int_sys.k = m.k;
    int_sys.s = con(i_con).s;
    int_sys.q = con(i_con).q;
    int_sys.h = con(i_con).h;
    
    int_sys.dydx = m.dydx;
    int_sys.dydu = m.dydu;
    
    int_sys.t = sol_sys.x;
    int_sys.x = @(t)devals(sol_sys, t);
    int_sys.u = con(i_con).u;
    int_sys.y = @(t)y(t, devals(sol_sys, t), u(t));
    
    int_sys.ie = sol_sys.ie;
    int_sys.te = sol_sys.xe;
    int_sys.xe = sol_sys.ye;
    int_sys.ue = u(int_sys.te);
    int_sys.ye = y(int_sys.te, int_sys.xe, int_sys.ue);
    
    int_sys.sol = sol_sys;
    
    int_sys.UseParams = opts.UseParams;
    int_sys.UseSeeds = opts.UseSeeds(:,i_con);
    
    % Distribute times for each observation
    int_sys = repmat(int_sys, n_obj,1);
    
    % Determine which objectives to evaluate for this experiment (those
    % that are not objectiveZero)
    isobjzero = strcmp('Objective.Data.Zero',{obj(:,i_con).Type});
    nonzeroobjs = find(~isobjzero);
    
    for i_obj = nonzeroobjs
        if obj(i_obj,i_con).Complex
            % Only reveal time points in range of observation
            % Note: deval will still not throw an error outside this range
            int_sys(i_obj).t = [int_sys(i_obj).t(int_sys(i_obj).t < obj(i_obj,i_con).tF), obj(i_obj,i_con).tF];
        else
            % Evaluate all requested time points
            int_sys(i_obj).t = obj(i_obj,i_con).DiscreteTimes;
            int_sys(i_obj).x = int_sys(i_obj).x(int_sys(i_obj).t);
            int_sys(i_obj).u = int_sys(i_obj).u(int_sys(i_obj).t);
            int_sys(i_obj).y = int_sys(i_obj).y(int_sys(i_obj).t);
        end
    end
    
    % *Compute G*
    % Extract continuous term
    if opts.continuous(i_con)
        G_cont = int_sys(1).sol.y(nx+1,end);
    else
        G_cont = 0;
    end
    
    % Compute discrete term
    G_disc = 0;
    discrete_times_all = cell(n_obj,1);
    for i_obj = nonzeroobjs
        [iDiscG, temp] = obj(i_obj,i_con).G(int_sys(i_obj));
        discrete_times_all{i_obj} = row(unique(temp));
        G_disc = G_disc + opts.ObjWeights(i_obj,i_con) * iDiscG;
    end
    
    discrete_times = vec(unique([discrete_times_all{:}]));
    
    % Add to cumulative goal value
    G = G + G_cont + G_disc;
    
    % * Integrate Adjoint *
    % Construct system
    [der, jac, del] = constructAdjointSystem();
    
    % Set initial conditions
    ic = zeros(nx+inT,1);
    
    % Integrate [lambda; D] backward in time
    sol = accumulateOdeRevSelect(der, jac, 0, tF, ic, [con(i_con).Discontinuities; discrete_times], 0, [], opts.RelTol, opts_i.AbsTol(nx+opts.continuous(i_con)+1:nx+opts.continuous(i_con)+nx+inT), del);
    
    % * Complete steady-state *
    if con(i_con).SteadyState
        % * Start Adjoint again *
        [der, jac] = constructSteadystateSystem();
        
        % Set initial conditions starting from end of previous run
        ic = sol.y;
        
        % Integrate [lambda; D] backward in time and replace previous run
        sol = accumulateOdeRevSelect(der, jac, 0, ssSol.xe, ic, con(i_con).private.BasalDiscontinuities, 0, [], opts.RelTol, opts_i.AbsTol(nx+opts.continuous(i_con)+1:nx+opts.continuous(i_con)+nx+inT));
    end
    
    % *Add contributions to derivative*
    % (Subtract contribution because gradient was integrated backward)
    
    % Rate parameters
    curD = zeros(nT,1);
    curD([1:nTk, Tsind+1:Tsind+inTs, Tqind+1:Tqind+inTq, Thind+1:Thind+inTh]) = -sol.y(nx+1:end,end);
    
    % Initial conditions
    lambda = -sol.y(1:nx,end);
    dx0ds_val = m.dx0ds(con(i_con).s);
    curD(Tsind+1:Tsind+inTs) = vec(curD(Tsind+1:Tsind+inTs)) + dx0ds_val(:,UseSeeds_i).' * lambda;
    
    % Add to cumulative goal value
    D = D + curD;
    
    % Update parameter index positions
    Tsind = Tsind + inTs;
    Tqind = Tqind + inTq;
    Thind = Thind + inTh;
    
    if verboseAll; fprintf('iCon = %d\t|dGdT| = %g\tTime = %0.2f\n', i_con, norm(curD), toc); end    
end

if opts.Verbose; fprintf('Summary: |dGdT| = %g\n', norm(D)); end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating lambda and D %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructAdjointSystem()
        dx0dd   = m.dx0ds;
        dfdx    = m.dfdx;
        dfdu    = m.dfdu;
        dfdk    = m.dfdk;
        dfdT    = @dfdTSub;
        dudq    = con(i_con).dudq;
        dddh    = con(i_con).dddh;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [lambda; D] with respect to time
        function val = derivative(t, joint)
            u_i = u(t);
            x_i = deval(sol_sys, t, 1:nx);
            y_i = y(t, x_i, u_i);
            dydx_i = dydx(t, x_i, u_i);
            l = joint(1:nx);
            
            % Sum continuous objective functions
            dgdx = zeros(nx,1);
            for i = 1:n_obj
                dgdx = dgdx + dydx_i.' * opts.ObjWeights(i,i_con)*obj(i,i_con).dgdy(t,y_i);
            end
            
            val = [dgdx; zeros(inT,1)] - [dfdx(t,x_i,u_i).'; dfdT(t,x_i,u_i).'] * l;
        end
        
        % Jacobian of [lambda; D] derivative
        function val = jacobian(t, joint)
            ui = u(t);
            x = deval(sol_sys, t, 1:nx);
            
            val = [-dfdx(t,x,ui).', sparse(nx,inT);
                   -dfdT(t,x,ui).', sparse(inT,inT)];
        end
        
        % Discrete effects of the objective function
        function val = delta(t, joint)
            u_i = u(t);
            x_i = deval(sol_sys, t, 1:nx);
            dydx_i = dydx(t, x_i, u_i);
            dydk_i = dydk(t, x_i, u_i);
            dGdx = zeros(nx,1);
            dGdT = zeros(inT,1);
            for i = 1:n_obj
                dGdy = obj(i,i_con).dGdy(t, int_sys(i));
                dGdx = dGdx + dydx_i.' * opts.ObjWeights(i,i_con)*dGdy;
                dGdk = obj(i,i_con).dGdk(int_sys(i)) + dydk_i.' * dGdy; % k_ % partial dGdk(i)
                dGds = obj(i,i_con).dGds(int_sys(i)); % s_ % partial dGds(i)
                dGdq = obj(i,i_con).dGdq(int_sys(i)); % q_ % partial dGdq(i)
                dGdh = obj(i,i_con).dGdh(int_sys(i)); % h_ % partial dGdh(i)
                dGdT = dGdT + opts.ObjWeights(i,i_con)*[dGdk(opts.UseParams); dGds(UseSeeds_i); dGdq(UseInputControls_i); dGdh(UseDoseControls_i)]; % T_ + (k_ -> T_) -> T_
            end
            
            lambda = -joint(1:nx,end) + dGdx; % Update current lambda
            dx0dd_i = dx0dd(d(t));
            dddh_i = dddh(t);
            dddh_i = dddh_i(:,UseDoseControls_i);
            dose_change = dddh_i.' * dx0dd_i.' * lambda;
            dGdT = dGdT + [zeros(nTk,1); zeros(inTs,1); zeros(inTq,1); dose_change];
            
            val = [dGdx; dGdT];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = dfdk(t,x,u);
            dfdq = dfdu(t,x,u) * dudq(t);
            val = [val(:,opts.UseParams), sparse(nx,inTs), dfdq(:,UseInputControls_i), sparse(nx,inTh)];
        end
    end

    function [der, jac] = constructSteadystateSystem()
        dfdx = m.dfdx;
        dfdu = m.dfdu;
        dfdk = m.dfdk;
        dfdT = @dfdTSub;
        basal_u = con(i_con).private.basal_u;
        basal_dudq = con(i_con).private.basal_dudq;
        
        der = @derivative;
        jac = @jacobian;
        
        % Derivative of [lambda; D] with respect to time
        function val = derivative(t, joint)
            ui = basal_u(t);
            x = deval(ssSol, t, 1:nx);
            l = joint(1:nx);
            
            val = -[dfdx(-1,x,ui).'; dfdT(t,x,ui).'] * l;
        end
        
        % Jacobian of [lambda; D] derivative
        function val = jacobian(t, joint)
            ui = basal_u(t);
            x = deval(ssSol, t, 1:nx);
            
            val = [-dfdx(-1,x,ui).', sparse(nx,inT);
                   -dfdT(t,x,ui).', sparse(inT,inT)];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = dfdk(-1,x,u);
            dfdq = dfdu(-1,x,u) * basal_dudq(t);
            val = [val(:,opts.UseParams), zeros(nx,inTs), dfdq(:,UseInputControls_i), sparse(nx,inTh)];
        end
    end

    function [der, jac, del] = constructObjectiveSystem()
        f       = m.f;
        dfdx    = m.dfdx;
        x0      = m.x0;
        nd      = m.ns;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [x; G] with respect to time
        function val = derivative(t, joint)
            ui = u(t);            
            x = joint(1:nx);
            
            % Sum continuous objective functions
            g = 0;
            for i = 1:n_obj
                g = g + opts.ObjWeights(i) * obj(i).g(t,x,ui);
            end
            
            val = [f(t,x,ui); g];
        end
        
        % Jacobian of [x; G] derivative
        function val = jacobian(t, joint)
            ui = u(t);
            x = joint(1:nx);
            
            % Sum continuous objective gradients
            dgdx = zeros(1,nx);
            for i = 1:n_obj
                dgdx = dgdx + opts.ObjWeights(i) * vec(obj(i).dgdx(t,x,ui)).';
            end
            
            val = [dfdx(t,x,ui), sparse(nx,1);
                          dgdx,            0];
        end

        % Dosing
        function val = delta(t, joint)
            val = [x0(d(t)) - x0(zeros(nd,1)); 0];
        end
    end

    function [der, jac, del] = constructSystem()
        f     = m.f;
        dfdx  = m.dfdx;
        x0    = m.x0;
        nd    = m.ns;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of x with respect to time
        function val = derivative(t, x)
            ui   = u(t);
            val = f(t,x,ui);
        end
        
        % Jacobian of x derivative
        function val = jacobian(t, x)
            ui   = u(t);
            val = dfdx(t,x,ui);
        end
        
        % Dosing
        function val = delta(t, x)
            val = x0(d(t)) - x0(zeros(nd,1));
        end
    end
end
