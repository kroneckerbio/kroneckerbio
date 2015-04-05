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
nCon = numel(con);
nObj = size(obj,1);

% Initialize variables
G = 0;
D = zeros(nT,1);
Tsind = nTk; % Stores the position in D where the first x0 parameter goes for each iCon
Tqind = nTk+nTs; % Stores the position in D where the first q parameter goes for each iCon
Thind = nTk+nTs+nTq; % Stores the position in D where the first h parameter goes for each iCon
intOpts = opts;

if opts.Verbose; disp('Integrating adjoint...'); end
for iCon = 1:nCon
    if verboseAll; tic; end
    
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.ObjWeights = opts.ObjWeights(:,iCon);

    UseSeeds_i = opts.UseSeeds(:,iCon);
    intOpts.UseSeeds = UseSeeds_i;
    inTs = nnz(UseSeeds_i);
    
    UseInputControls_i = opts.UseInputControls{iCon};
    intOpts.UseControls = UseInputControls_i;
    inTq = nnz(UseInputControls_i);
    
    UseDoseControls_i = opts.UseDoseControls{iCon};
    intOpts.UseControls = UseDoseControls_i;
    inTh = nnz(UseDoseControls_i);
    
    inT = nTk + inTs + inTq + inTh;
    
    % Seeds
    s = con(iCon).s;
    
    % Input
    u_f = con(iCon).u;
    d_f = con(iCon).d;
    q = con(iCon).q;
    h = con(iCon).h;
    
    % * Integrate to steady-state
    if con(iCon).SteadyState
        ssSol = integrateSteadystateSys(m, con(iCon), intOpts);
        
        % Apply steady-state solution to initial conditions
        ic = ssSol.y(:,end);
    else
        ic = m.dx0ds * s + m.x0c;
    end
    
    % * Integrate system *
    % TODO: consider reusing the existing methods for doing this
    % Do not use select methods since the solution is needed at all time
    if opts.continuous(iCon)
        [der, jac, del] = constructObjectiveSystem();
        xSol = accumulateOdeFwd(der, jac, 0, con(iCon).tF, [ic; 0], con(iCon).Discontinuities, 1:nx, opts.RelTol, opts.AbsTol{iCon}(1:nx+1), del);
    else
        [der, jac, del] = constructSystem();
        xSol = accumulateOdeFwd(der,jac, 0, con(iCon).tF, ic, con(iCon).Discontinuities, 1:nx, opts.RelTol, opts.AbsTol{iCon}(1:nx), del);
    end
    xSol.nx = nx;
    xSol.u = con(iCon).u;
    xSol.y_ = m.y;
    xSol.dydx = m.dydx;
    xSol.dydu = m.dydu;
    xSol.k = m.k;
    xSol.s = s;
    xSol.q = q;
    xSol.h = h;
    xSol.UseParams = opts.UseParams;
    xSol.UseSeeds = UseSeeds_i;
    xSol.UseInputControls = UseInputControls_i;
    xSol.UseDoseControls = UseDoseControls_i;
    
    % Extract continuous term
    if opts.continuous(iCon)
        contG = xSol.y(nx+1,end);
    else
        contG = 0;
    end
    
    % Compute discrete term
    discG = 0;
    discreteTimes = [];
    for iObj = 1:nObj
        [iDiscG, temp] = obj(iObj,iCon).G(xSol);
        discreteTimes = [discreteTimes; vec(temp)];
        discG = discG + opts.ObjWeights(iObj,iCon) * iDiscG;
    end
    
    % Remove repetitive discreteTimes
    discreteTimes = unique(discreteTimes);
    nDisc = numel(discreteTimes);

    % Add to cumulative goal value
    G = G + contG + discG;
    
    % * Integrate Adjoint *
    % Construct system
    [der, jac, del] = constructAdjointSystem();
    
    % Set initial conditions
    ic = zeros(nx+inT,1);
    
    % Integrate [lambda; D] backward in time
    sol = accumulateOdeRevSelect(der, jac, 0, con(iCon).tF, ic, [con(iCon).Discontinuities; discreteTimes], 0, [], opts.RelTol, opts.AbsTol{iCon}(nx+opts.continuous(iCon)+1:nx+opts.continuous(iCon)+nx+nT), del);
    
    % * Complete steady-state *
    if con(iCon).SteadyState
        % * Start Adjoint again *
        [der, jac] = constructSteadystateSystem();
        
        % Set initial conditions starting from end of previous run
        ic = sol.y;
        
        % Integrate [lambda; D] backward in time and replace previous run
        sol = accumulateOdeRevSelect(der, jac, 0, ssSol.x(end), ic, [], 0, [], opts.RelTol, opts.AbsTol{iCon}(nx+opts.continuous(iCon)+1:nx+opts.continuous(iCon)+nx+nT));
    end
    
    % *Add contributions to derivative*
    % (Subtract contribution because gradient was integrated backward)
    
    % Rate parameters
    curD = zeros(nT,1);
    curD([1:nTk, Tsind+1:Tsind+inTs, Tqind+1:Tqind+inTq, Thind+1:Thind+inTh]) = -sol.y(nx+1:end,end);
    
    % Initial conditions
    lambda = -sol.y(1:nx,end);
    curD(Tsind+1:Tsind+inTs) = curD(Tsind+1:Tsind+inTs) + m.dx0ds(:,UseSeeds_i).' * lambda;
    
    % Add to cumulative goal value
    D = D + curD;
    
    % Update parameter index positions
    Tsind = Tsind + inTs;
    Tqind = Tqind + inTq;
    Thind = Thind + inTh;
    
    if verboseAll; fprintf('iCon = %d\t|dGdT| = %g\tTime = %0.2f\n', iCon, norm(curD), toc); end    
end

if opts.Verbose; fprintf('Summary: |dGdT| = %g\n', norm(D)); end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating lambda and D %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructAdjointSystem()
        dx0ds = m.dx0ds;
        dfdx = m.dfdx;
        dfdu = m.dfdu;
        dfdk = m.dfdk;
        dfdT = @dfdTSub;
        dudq = con(iCon).dudq;
        dddh = con(iCon).dddh;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [lambda; D] with respect to time
        function val = derivative(t, joint)
            u = u_f(t);
            x = deval(xSol, t, 1:nx);
            l = joint(1:nx);
            
            % Sum continuous objective functions
            dgdx = zeros(nx,1);
            for iObj = 1:nObj
                dgdx = dgdx + opts.ObjWeights(iObj,iCon)*obj(iObj,iCon).dgdx(t,x,u);
            end
            
            val = [dgdx; zeros(inT,1)] - [dfdx(t,x,u).'; dfdT(t,x,u).'] * l;
        end
        
        % Jacobian of [lambda; D] derivative
        function val = jacobian(t, joint)
            u = u_f(t);
            x = deval(xSol, t, 1:nx);
            
            val = [-dfdx(t,x,u).', sparse(nx,inT);
                   -dfdT(t,x,u).', sparse(inT,inT)];
        end
        
        % Discrete effects of the objective function
        function val = delta(t, joint)
            dGdx = zeros(nx,1);
            dGdT = zeros(inT,1);
            for iObj = 1:nObj
                dGdx = dGdx + opts.ObjWeights(iObj,iCon)*obj(iObj,iCon).dGdx(t, xSol);
                dGdk = obj(iObj,iCon).dGdk(t, xSol); % k_ % partial dGdk(i)
                dGds = obj(iObj,iCon).dGds(t, xSol); % s_ % partial dGds(i)
                dGdq = obj(iObj,iCon).dGdq(t, xSol); % q_ % partial dGdq(i)
                dGdh = obj(iObj,iCon).dGdh(t, xSol); % h_ % partial dGdq(i)
                dGdT = dGdT + opts.ObjWeights(iObj,iCon)*[dGdk(opts.UseParams); dGds(UseSeeds_i); dGdq(UseInputControls_i); dGdh(UseDoseControls_i)]; % T_ + (k_ -> T_) -> T_
            end
            
            lambda = -joint(1:nx,end) + dGdx; % Update current lambda
            dddh_i = dddh(t);
            dddh_i = dddh_i(:,UseDoseControls_i);
            dose_change = dddh_i.' * dx0ds.' * lambda;
            dGdT = dGdT + [zeros(nTk,1); zeros(nTs,1); zeros(nTq,1); dose_change];
            
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
        if opts.UseModelInputs
            dudq = m.dudq;
        else
            dudq = con.dudq;
        end
        
        der = @derivative;
        jac = @jacobian;
        
        % Derivative of [lambda; D] with respect to time
        function val = derivative(t, joint)
            u = u_f(-1);
            x = deval(ssSol, t, 1:nx);
            l = joint(1:nx);
            
            val = -[dfdx(-1,x,u).'; dfdT(-1,x,u).'] * l;
        end
        
        % Jacobian of [lambda; D] derivative
        function val = jacobian(t, joint)
            u = u_f(-1);
            x = deval(ssSol, t, 1:nx);
            
            val = [-dfdx(-1,x,u).', sparse(nx,inT);
                   -dfdT(-1,x,u).', sparse(inT,inT)];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = dfdk(-1,x,u);
            dfdq = dfdu(-1,x,u) * dudq(-1);
            val = [val(:,opts.UseParams), zeros(nx,inTs), dfdq(:,UseInputControls_i), sparse(nx,inTh)];
        end
    end

    function [der, jac, del] = constructObjectiveSystem()
        f     = m.f;
        dfdx  = m.dfdx;
        dx0ds = m.dx0ds;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [x; G] with respect to time
        function val = derivative(t, joint)
            u = u_f(t);            
            x = joint(1:nx);
            
            % Sum continuous objective functions
            g = 0;
            for i = 1:nObj
                g = g + opts.ObjWeights(i) * obj(i).g(t,x,u);
            end
            
            val = [f(t,x,u); g];
        end
        
        % Jacobian of [x; G] derivative
        function val = jacobian(t, joint)
            u = u_f(t);
            x = joint(1:nx);
            
            % Sum continuous objective gradients
            dgdx = zeros(1,nx);
            for i = 1:nObj
                dgdx = dgdx + opts.ObjWeights(i) * vec(obj(i).dgdx(t,x,u)).';
            end
            
            val = [dfdx(t,x,u), sparse(nx,1);
                          dgdx,            0];
        end

        % Dosing
        function val = delta(t, joint)
            val = [dx0ds * d_f(t); 0];
        end
    end

    function [der, jac, del] = constructSystem()
        f     = m.f;
        dfdx  = m.dfdx;
        dx0ds = m.dx0ds;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of x with respect to time
        function val = derivative(t, x)
            u   = u_f(t);
            val = f(t,x,u);
        end
        
        % Jacobian of x derivative
        function val = jacobian(t, x)
            u   = u_f(t);
            val = dfdx(t,x,u);
        end
        
        % Dosing
        function val = delta(t, x)
            val = dx0ds * d_f(t);
        end
    end
end
