function [G, D] = computeObjSensAdj(m, con, obj, opts)
verboseAll = max(opts.Verbose-1,0);
if verboseAll; tic; end

% Constants
nx = m.nx;
nTk = sum(opts.UseParams);
nTx = sum(sum(opts.UseICs));
nTq = sum(cat(1,opts.UseControls{:}));
nT  = nTk + nTx + nTq;
nCon = numel(con);
nObj = size(obj,1);

% Initialize variables
G = 0;
D = zeros(nT,1);
Txind = nTk; % Stores the position in D where the first x0 parameter goes for each iCon
Tqind = nTk+nTx; % Stores the position in D where the first q parameter goes for each iCon
intOpts = opts;

if opts.Verbose; disp('Integrating adjoint...'); end
for iCon = 1:nCon
    if verboseAll; tic; end
    
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.ObjWeights = opts.ObjWeights(:,iCon);

    % If opts.UseModelICs is false, the number of variables can change
    if opts.UseModelICs
        inTx = nTx;
    else
        intOpts.UseICs = opts.UseICs(:,iCon);
        inTx = sum(intOpts.UseICs);
    end
    
    % If opts.UseModelInputs is false, the number of variables can change
    if opts.UseModelInputs
        inTq = nTq;
    else
        intOpts.UseControls = opts.UseControls(iCon);
        inTq = sum(intOpts.UseControls{1});
    end
    
    inT = nTk + inTx + inTq;
    
    % * Integrate to steady-state
    if con(iCon).SteadyState
        ss = true;
        ssSol = integrateSteadystateSys(m, con(iCon), intOpts);
        
        % Apply steady-state solution to initial conditions
        if opts.UseModelICs
            m = m.Update(m.k, ssSol.y(:,end), m.q);
        else
            con(iCon) = con(iCon).Update(ssSol.y(:,end), con(iCon).q);
        end
        
        % Don't do it again
        con(iCon).SteadyState = false;
    else
        ss = false;
    end
    
    % * Integrate system *
    % Do not use select methods since the solution is needed at all time
    if opts.continuous(iCon)
        xSol = integrateObj(m, con(iCon), obj(:,iCon), intOpts);
    else
        xSol = integrateSys(m, con(iCon), intOpts);
    end
    
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
    [der, jac, del] = constructSystem();
    
    % Set initial conditions
    ic = zeros(nx+nTk+inTq,1);
    
    % Input
    if opts.UseModelInputs
        u = m.u;
    else
        u = con(iCon).u;
    end
    
    % Integrate [lambda; D] backward in time
    sol = accumulateOde(der, jac, 0, con(iCon).tF, ic, u, [con(iCon).Discontinuities; discreteTimes], [], opts.RelTol, opts.AbsTol{iCon}(nx+opts.continuous(iCon)+1:nx+opts.continuous(iCon)+nx+nTk+inTq), del, -1, [], [], [], 0);
    
    % * Complete steady-state *
    if ss
        % * Start Adjoint again *
        [der, jac] = constructSteadystateSystem();
        
        % Set initial conditions starting from end of previous run
        ic = sol.y;
        
        % Integrate [lambda; D] backward in time and replace previous run
        sol = accumulateOde(der, jac, 0, ssSol.x(end), ic, u, [], [], opts.RelTol, opts.AbsTol{iCon}(nx+opts.continuous(iCon)+1:nx+opts.continuous(iCon)+nx+nTk+inTq), [], -1, [], [], [], 0);
    end
    
    % *Add contributions to derivative*
    % (Subtract contribution because gradient was integrated backward)
    
    % Rate parameters
    curD = zeros(nT,1);
    curD([1:nTk, Tqind+1:Tqind+inTq]) = -sol.y(nx+1:end,end);
    
    % Initial conditions
    if opts.UseModelICs
        lambda = sol.y(1:nx,end);
        dx0dx0 = speye(nx,nx);
        curD(Txind+1:Txind+inTx) = dx0dx0(opts.UseICs,:) * -lambda;
    else
        lambda = sol.y(1:nx,end);
        dx0dx0 = speye(nx,nx);
        curD(Txind+1:Txind+inTx) = dx0dx0(opts.UseICs(:,iCon),:) * -lambda;
    end
    
    % Add to cumulative goal value
    D = D + curD;
    
    % Update condition x0 position
    if ~opts.UseModelICs
        Txind = Txind + inTx;
    end
    % Update condition q position
    if ~opts.UseModelInputs
        Tqind = Tqind + inTq;
    end
    
    if verboseAll; fprintf('iCon = %d\t|dGdT| = %g\tTime = %0.2f\n', iCon, norm(curD), toc); end    
end

if opts.Verbose; fprintf('Summary: |dGdT| = %g\n', norm(D)); end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating lambda and D %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        dfdx = m.dfdx;
        dfdu = m.dfdu;
        dfdT = @dfdTSub;
        if opts.UseModelInputs
            dudq = m.dudq;
        else
            dudq = con.dudq;
        end
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [lambda; D] with respect to time
        function val = derivative(t, joint, u)
            u = u(t);
            x = deval(xSol, t, 1:nx);
            l = joint(1:nx);
            
            % Sum continuous objective functions
            dgdx = zeros(nx,1);
            dgdT = zeros(nTk+inTq,1);
            for iObj = 1:nObj
                dgdx = dgdx + opts.ObjWeights(iObj,iCon)*vec(obj(iObj,iCon).dgdx(t,x,u));
                dgdk = obj(iObj,iCon).dgdk(t,x,u); % k_
                dgdT = dgdT + opts.ObjWeights(iObj,iCon)*[vec(dgdk(opts.UseParams)); zeros(inTq,1)]; % T_ + (k_ -> T_) -> T_
            end
            
            val = [dgdx; dgdT] - [dfdx(t,x,u).'; dfdT(t,x,u).'] * l;
        end
        
        % Jacobian of [lambda; D] derivative
        function val = jacobian(t, joint, u)
            u = u(t);
            x = deval(xSol, t, 1:nx);
            
            val = [-dfdx(t,x,u).', sparse(nx,nTk+inTq);
                   -dfdT(t,x,u).', sparse(nTk+inTq,nTk+inTq)];
        end
        
        % Discrete effects of the objective function
        function val = delta(t)
            dGdx = zeros(nx,1);
            dGdT = zeros(nTk+inTq,1);
            for iObj = 1:nObj
                dGdx = dGdx + opts.ObjWeights(iObj,iCon)*vec(obj(iObj,iCon).dGdx(t, xSol));
                dGdk = obj(iObj,iCon).dGdk(t, xSol); % k_
                dGdT = dGdT + opts.ObjWeights(iObj,iCon)*[vec(dGdk(opts.UseParams)); zeros(inTq,1)]; % T_ + (k_ -> T_) -> T_
            end
            
            val = [dGdx; dGdT];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = m.dfdk(t,x,u);
            dfdq = dfdu(t,x,u) * dudq(t);
            val = [val(:,opts.UseParams) dfdq(:,opts.UseControls{1})];
        end
    end

    function [der, jac] = constructSteadystateSystem()
        dfdx = m.dfdx;
        dfdu = m.dfdu;
        dfdT = @dfdTSub;
        if opts.UseModelInputs
            dudq = m.dudq;
        else
            dudq = con.dudq;
        end
        
        der = @derivative;
        jac = @jacobian;
        
        % Derivative of [lambda; D] with respect to time
        function val = derivative(t, joint, u)
            u = u(-1);
            x = deval(ssSol, t, 1:nx);
            l = joint(1:nx);
            
            val = -[dfdx(-1,x,u).'; dfdT(-1,x,u).'] * l;
        end
        
        % Jacobian of [lambda; D] derivative
        function val = jacobian(t, joint, u)
            u = u(-1);
            x = deval(ssSol, t, 1:nx);
            
            val = [-dfdx(-1,x,u).', sparse(nx,nTk+inTq);
                   -dfdT(-1,x,u).', sparse(nT,nTk+inTq)];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = m.dfdk(-1,x,u);
            dfdq = dfdu(-1,x,u) * dudq(-1);
            val = [val(:,opts.UseParams) dfdq(:,opts.UseControls{1})];
        end
    end
end