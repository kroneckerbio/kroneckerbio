function sol = integrateDblsensSelect(m, con, tGet, opts)
% Constants
nx = m.nx;
nu = m.nu;
nq = numel(opts.UseControls{1});
nk = m.nk;
nTk = sum(opts.UseParams);
nTx = sum(sum(opts.UseICs));
nTq = sum(opts.UseControls{1});
nT  = nTk + nTx + nTq;

% Construct system
[der, jac] = constructSystem();

if ~con.SteadyState
    % Initial conditions [x0; vec(dxdT0)]
    if opts.UseModelICs
        x0 = m.x0;
    else
        x0 = con.x0;
    end
    
    % Initial effect of rates on sensitivities is 0
    dxdT0 = zeros(nx, nTk); % Active rate parameters
    
    % Initial effect of ics on sensitivities is 1 for that state
    dxdTx0                = zeros(nx,nTx);
    dxdTx0(opts.UseICs,:) = eye(nTx);
    
    % Initial effect of qs on sensitivities is 0
    dxdTq = zeros(nx, nTq);
    
    % Initial double sensitivities are zero
    d2xdT2 = zeros(nx*nT*nT,1);
    
    % Combine them into a vector
    ic = [x0; vec([dxdT0, dxdTx0, dxdTq]); d2xdT2];
else
    % Run to steady-state first
    ic = steadystateDblsens(m, con, opts);
end

% Input
if opts.UseModelInputs
    u = m.u;
else
    u = con.u;
end

% Integrate [x; dxdT; d2xdT2] with respect to time
sol = accumulateOde(der, jac, 0, con.tF, ic, u, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+nx*nT+nx*nT*nT), [], 1, [], [], [], tGet);
sol.u = u(tGet);
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x and dxdT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac] = constructSystem()
        
        IT = speye(nT);
        IT2 = speye(nT*nT);
        TPermuteInd = vec(permute(reshape(1:nx*nT*nT, nx,nT,nT), [1,3,2])); % indexes that permute the last two T of fTT_
        fkUseParams = reshape(1:nx*nTk, nx,nTk);
        fkUseParams = vec(fkUseParams(:,opts.UseParams)); % indexes that remove inactive k from fk_
%         fqUseControls = reshape(1:nx*nTq, nx,nTq);
%         fqUseControls = vec(fqUseControls(:,opts.UseControls{1})); % indexes that remove inactive q from fq_

        dxdTStart = nx+1;
        dxdTEnd   = nx+nx*nT;
        d2xdT2Start = nx+nx*nT+1;
        d2xdT2End   = nx+nx*nT+nx*nT*nT;
        
        f       = m.f;
        dfdx    = m.dfdx;
        dfdu    = m.dfdu;
        dfdT    = @dfdTSub;
        d2fdx2  = m.d2fdx2;
        d2fdT2  = @d2fdT2Sub;
        d2fdudx = m.d2fdudx;
        d2fdTdx = @d2fdTdxSub;
        d2fdxdT = @d2fdxdTSub;
        if opts.UseModelInputs
            dudq = m.dudq;
        else
            dudq = con.dudq;
        end
        
        der = @derivative;
        jac = @jacobian;
        
        % Derivative of [x; dxdT; d2xdT2] with respect to time
        function val = derivative(t, joint, u)
            u = u(t);            
            x = joint(1:nx);
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % xT_ -> x_T
            d2xdT2 = reshape(joint(d2xdT2Start:d2xdT2End), nx,nT*nT); % xTT_ -> x_TT
            
            % Derivative of dxdT
            % dxdTdot = dfdx *{x.x} dxdT + dxdT
            dxdTdot = dfdx(t,x,u) * dxdT + dfdT(t,x,u); % f_x * x_T + f_T -> f_T
            
            % Derivative of d2xdT2
            % d2xdT2dot = dfdx *{x.x} d2xdT1dT2 + (d2fdx2 *{x.x} dxdT2) *{x.x} dxdT1 + d2fdT2dx *{x.x} dxdT1 + d2fdT1dx *{x.x} dxdT2 + d2fdT1dT2
            temp = sparse(d2fdxdT(t,x,u) * dxdT); % fT_x * x_T -> fT_T
            temp = temp + spermute132(temp, [nx,nT,nT], [nx*nT,nT]); % fT_T + (fT_T -> fT_T) -> fT_T
            d2xdT2dot = sparse(d2fdx2(t,x,u) * dxdT); % fx_x * x_T -> fx_T
            d2xdT2dot = spermute132(d2xdT2dot, [nx,nx,nT], [nx*nT,nx]); % fx_T -> fT_x
            d2xdT2dot = sparse(d2xdT2dot * dxdT); %fT_x * x_T -> fT_T
            d2xdT2dot = reshape(dfdx(t,x,u) * d2xdT2, nx*nT,nT) + d2xdT2dot + temp + d2fdT2(t,x,u); % (f_x * x_TT -> f_TT -> fT_T) + fT_T + fT_T + fT_T -> fT_T
            
            % Combine
            val = [f(t,x,u); vec(dxdTdot); vec(d2xdT2dot)];
        end
        
        % Jacobian of [x; dxdT; d2xdT2] derivative
        function val = jacobian(t, joint, u)
            u = u(t);            
            x = joint(1:nx); % x_
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % x_T
            d2xdT2 = reshape(joint(d2xdT2Start:d2xdT2End), nx,nT*nT); % xTT_ -> x_TT
            
            % Compute d/dx(dfdT)
            d2fdxdTVal = sparse(d2fdx2(t,x,u) * dxdT) + d2fdTdx(t,x,u); % fx_T
            d2fdxdTVal = spermute132(d2fdxdTVal, [nx,nx,nT], [nx*nT,nx]);
            
            % Compute d/dx(d2fdT2)
            d3fdxdT2 = sparse(nx*nT*nT,nx);% TEMP
            
            % Compute d/dxdT(d2fdT2)
            % d/dxdT(d2fdT2) = (d2fdx2 *{x.x} dxdT2) * IT1 + (d2fdx2 *{x.x} dxdT1) * IT2 + d2fdT2dx * IT1 + d2fdT1dx * IT2
            d3fddxdTT2 = sparse(d2fdx2(t,x,u) * dxdT); % fx_x * x_T -> fx_T
            d3fddxdTT2 = spermute132(d3fddxdTT2, [nx,nx,nT], [nx*nT,nx]); % fx_T -> fT_x
            d3fddxdTT2 = d3fddxdTT2 + d2fdxdT(t,x,u); % fT_x + fT_x -> fT_x
            d3fddxdTT2 = kron(IT, d3fddxdTT2); % T_T kron fT_x -> fTT_xT
            d3fddxdTT2 = d3fddxdTT2 + d3fddxdTT2(TPermuteInd,:); % fTT_xT + (fTT_xT -> fTT_xT) -> fTT_xT
            
            % Combine
            val = [dfdx(t,x,u), sparse(nx,nx*nT),      sparse(nx,nx*nT*nT);
                   d2fdxdTVal,  kron(IT, dfdx(t,x,u)), sparse(nx*nT,nx*nT*nT);
                   d3fdxdT2,    d3fddxdTT2,            kron(IT2, dfdx(t,x,u))];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = m.dfdk(t,x,u);
            dfdq = dfdu(t,x,u) * dudq(t);
            val = [val(:,opts.UseParams) sparse(nx, nTx) dfdq(:,opts.UseControls{1})];
        end
        
        % Modifies d2fdkdk to relate only to the parameters of interest
        function val = d2fdT2Sub(t, x, u)
            val = m.d2fdk2(t,x,u);
            
            d2fdq2 = m.d2fdu2(t,x,u) * dudq(t); % fu_u * u_q -> fu_q
            d2fdq2 = d2fdq2(:,opts.UseControls{1}); % fu_q -> fu_q(T)
            d2fdq2 = spermute132(d2fdq2, [nx,nu,nTq], [nx*nTq,nu]); % (fu_q -> fq_u) * u_q -> fq_q
            d2fdq2 = d2fdq2(:,opts.UseControls{1}); % fq_q -> fq_q(T)
            
            d2fdkdq = m.d2fdudk(t,x,u) * dudq(t); % fk_u * u_q -> fk_q
            d2fdkdq = d2fdkdq(:,opts.UseControls{1}); % fk_q -> fk_q(T)
            d2fdkdq = spermute132(d2fdkdq, [nx,nk,nTq], [nx*nTq,nk]); % fk_q -> fq_k
            d2fdkdq = d2fdkdq(:,opts.UseParams); % fq_k -> fq_k(T)
            
            val = [val(fkUseParams,opts.UseParams), sparse(nx*nTk, nTx+nTq), spermute132(d2fdkdq, [nx,nTq,nTk], [nx*nTk,nTq]);
                   sparse(nx*nTx, nTk+nTx+nTq);
                   d2fdkdq,                         sparse(nx*nTq, nTx),     d2fdq2];
        end
        
        % Modifies d2fdkdx to relate only to the parameters of interest
        function val = d2fdTdxSub(t, x, u)
            val = m.d2fdkdx(t,x,u); % fx_k
            d2fdqdx = d2fdudx(t,x,u) * dudq(t); % fx_u * u_q -> fx_q
            val = [val(:,opts.UseParams) sparse(nx*nx, nTx) d2fdqdx(:,opts.UseControls{1})]; % fx_T
        end
        
        % Modifies d2fdxdk to relate only to the parameters of interest
        function val = d2fdxdTSub(t, x, u)
            val = spermute132(d2fdTdx(t,x,u), [nx,nx,nT], [nx*nT,nx]); % fx_T -> fT_x
        end
end

end