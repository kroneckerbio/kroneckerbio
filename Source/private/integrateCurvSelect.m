function sol = integrateCurvSelect(m, con, tGet, opts)
% Constants
nx = m.nx;
nu = m.nu;
nk = m.nk;
ns = m.ns;
nq = con.nq;
nh = con.nh;
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
    dx0dTk = zeros(nx, nTk); % Active rate parameters
    
    % Initial effect of seeds on states is dx0ds
    dx0dTs = m.dx0ds(:,opts.UseSeeds);
    
    % Initial effect of qs on sensitivities is 0
    dx0dTq = zeros(nx, nTq);
    
    % Initial effect of hs on sensitivities is 0
    dx0dTh = zeros(nx, nTh);

    % Initial curavtures are zero
    d2x0dT2 = zeros(nx*nT*nT,1);
    
    % Combine them into a vector
    ic = [x0; vec([dx0dTk, dx0dTs, dx0dTq, dx0dTh]); d2x0dT2];
else
    % Run to steady-state first
    ic = steadystateCurv(m, con, opts);
end

% Integrate [f; dfdT; d2fdT2] over time
sol = accumulateOdeFwdSelect(der, jac, 0, con.tF, ic, con.Discontinuities, tGet, 1:nx, opts.RelTol, opts.AbsTol(1:nx+nx*nT+nx*nT*nT), del);
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f and dfdT and d2fdT2 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        
        IT = speye(nT);
        IT2 = speye(nT*nT);
        TPermuteInd = vec(permute(reshape(1:nx*nT*nT, nx,nT,nT), [1,3,2])); % indexes that permute the last two T of fTT_
%         fkUseParams = reshape(1:nx*nk, nx,nk);
%         fkUseParams = vec(fkUseParams(:,opts.UseParams)); % indexes that remove inactive k from fk_
        fkUseParams = linearslicer([nx,nk], true(nx,1), opts.UseParams);
        fqUseInputControls = linearslicer([nx,nq], true(nx,1), opts.UseInputControls);
        uqUseInputControls = linearslicer([nu,nq], true(nu,1), opts.UseInputControls);
        dhUseDoseControls = linearslicer([ns,nh], true(ns,1), opts.UseDoseControls);

        dxdTStart = nx+1;
        dxdTEnd   = nx+nx*nT;
        d2xdT2Start = nx+nx*nT+1;
        d2xdT2End   = nx+nx*nT+nx*nT*nT;
        
        f       = m.f;
        dfdx    = m.dfdx;
        dfdu    = m.dfdu;
        dfdk    = m.dfdk;
        dfdT    = @dfdTSub;
        d2fdx2  = m.d2fdx2;
        d2fdT2  = @d2fdT2Sub;
        d2fdudx = m.d2fdudx;
        d2fdTdx = @d2fdTdxSub;
        d2fdxdT = @d2fdxdTSub;
        d2fdu2  = m.d2fdu2;
        d2fdk2  = m.d2fdk2;
        d2fdudk = m.d2fdudk;
        d2fdkdx = m.d2fdkdx;
        uf      = con.u;
        dudq    = con.dudq;
        d2udq2  = con.d2udq2;
        d       = con.d;
        dddh    = con.dddh;
        d2ddh2  = con.d2ddh2;
        dx0ds   = m.dx0ds;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [x; dxdT; d2xdT2] with respect to time
        function val = derivative(t, joint)
            u = uf(t);            
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
            d2xdT2dot = sparse(d2xdT2dot * dxdT); % fT_x * x_T -> fT_T
            d2xdT2dot = reshape(dfdx(t,x,u) * d2xdT2, nx*nT,nT) + d2xdT2dot + temp + d2fdT2(t,x,u); % (f_x * x_TT -> f_TT -> fT_T) + fT_T + fT_T + fT_T -> fT_T
            
            % Combine
            val = [f(t,x,u); vec(dxdTdot); vec(d2xdT2dot)];
        end
        
        % Jacobian of [x; dxdT; d2xdT2] derivative
        function val = jacobian(t, joint)
            u = uf(t);            
            x = joint(1:nx); % x_
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % x_T
            d2xdT2 = reshape(joint(d2xdT2Start:d2xdT2End), nx,nT*nT); % xTT_ -> x_TT
            
            % Compute d/dx(dfdT)
            d2fdxdT_val = sparse(d2fdx2(t,x,u) * dxdT) + d2fdTdx(t,x,u); % fx_T
            d2fdxdT_val = spermute132(d2fdxdT_val, [nx,nx,nT], [nx*nT,nx]);
            
            % Compute d/dx(d2fdT2)
            % This part is very computationally expensive and for unknown
            % reasons does not matter at all. Therefore, it is set to zero.
            % The computational expense mostly comes from sparse*full
            % yielding full because Matlab does not understand how sparse
            % these things really are.
            d3fdxdTdT_val = sparse(nx*nT*nT,nx);
            
            % Real computation
            %temp = sparse(d3fdxdTdx(t,x,u) * dxdT); % fxT_x * x_T -> fxT_T
            %temp = temp + spermute132(temp, [nx*nx,nT,nT], [nx*nx*nT,nT]); % fxT_T + (fxT_T -> fxT_T) -> fxT_T
            %d3fdxdTdT_val = sparse(d3fdx3(t,x,u) * dxdT); % fxx_x * x_T -> fxx_T
            %d3fdxdTdT_val = spermute132(d3fdxdTdT_val, [nx*nx,nx,nT], [nx*nx*nT,nx]); % fxx_T -> fxT_x
            %d3fdxdTdT_val = sparse(d3fdxdTdT_val * dxdT); % fxT_x * x_T -> fxT_T
            %d3fdxdTdT_val = reshape(df2dx2(t,x,u) * d2xdT2, nx*nx*nT,nT) + d3fdxdTdT_val + temp + d3fdxdTdT; % ((fx_x * x_TT -> fx_TT) -> fxT_T) + fxT_T + fxT_T + fxT_T -> fxT_T
            
            % Compute d/dxdT(d2fdT2)
            % d/dxdT(d2fdT2) = (d2fdx2 *{x.x} dxdT2) * IT1 + (d2fdx2 *{x.x} dxdT1) * IT2 + d2fdT2dx * IT1 + d2fdT1dx * IT2
            d3fddxdTT2 = sparse(d2fdx2(t,x,u) * dxdT); % fx_x * x_T -> fx_T
            d3fddxdTT2 = spermute132(d3fddxdTT2, [nx,nx,nT], [nx*nT,nx]); % fx_T -> fT_x
            d3fddxdTT2 = d3fddxdTT2 + d2fdxdT(t,x,u); % fT_x + fT_x -> fT_x
            d3fddxdTT2 = kron(IT, d3fddxdTT2); % T_T kron fT_x -> fTT_xT
            d3fddxdTT2 = d3fddxdTT2 + d3fddxdTT2(TPermuteInd,:); % fTT_xT + (fTT_xT -> fTT_xT) -> fTT_xT
            
            % Combine
            val = [dfdx(t,x,u),   sparse(nx,nx*nT),      sparse(nx,nx*nT*nT);
                   d2fdxdT_val,   kron(IT, dfdx(t,x,u)), sparse(nx*nT,nx*nT*nT);
                   d3fdxdTdT_val, d3fddxdTT2,            kron(IT2, dfdx(t,x,u))];
        end
        
        % Dosing
        function val = delta(t, joint)
            deltax = dx0ds * d(t);
            
            % dxdh = dxdd *{d.d} dddh
            dddh_i = dddh(t); % s_h
            dxdTh = dx0ds * dddh_i(:,opts.UseDoseControls); % x_s * (s_h -> s_H) -> x_H
            dxdT = [zeros(nx,nTk+nTs+nTq), dxdTh];
            
            % d2xdh2 = (dxdd2dd1 *{d.d} dd2dh2) *{d.d} dddh1 + dxdd *{d.d} d2ddh2dh1
            % Currently, first term is gauranteed to be zero because x0 is linear
            d2dh2_i = d2ddh2(t); %sh_h
            d2xTh2 = dx0ds * reshape(d2dh2_i(dhUseDoseControls, opts.UseDoseControls), ns,nTh*nTh); % x_s * (sh_h -> sH_H -> s_HH) -> x_HH
            d2xdT2 = [sparse(nx*(nTk+nTs+nTq),nT); sparse(nx*nTh,nTk+nTs+nTq), reshape(d2xTh2, [nx*nTh,nTh])];
            
            val = [deltax; vec(dxdT); vec(d2xdT2)];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = dfdk(t,x,u); % f_k
            dudq_i = dudq(t); % u_q
            dfdTq = dfdu(t,x,u) * dudq_i(:,opts.UseInputControls); % f_u * (u_q -> u_Q) -> f_Q
            val = [val(:,opts.UseParams), sparse(nx,nTs), dfdTq, sparse(nx,nTh)]; % f_T
        end
        
        % Modifies d2fdkdk to relate only to the parameters of interest
        function val = d2fdT2Sub(t, x, u)
            val = d2fdk2(t,x,u); % fk_k
            
            dudq_i = dudq(t); % u_q
            dudq_i = dudq_i(:,opts.UseDoseControls); % u_q -> u_Q
            
            d2udQ2_i = d2udq2(t); % uq_q
            d2udQ2_i = reshape(d2udQ2_i(uqUseInputControls,opts.UseInputControls), [nu,nTq*nTq]); % uq_q -> uQ_Q -> u_QQ
            
            d2fdq2 = d2fdu2(t,x,u) * dudq_i; % fu_u * u_Q -> fu_Q
            d2fdq2 = spermute132(d2fdq2, [nx,nu,nTq], [nx*nTq,nu]) * dudq_i + reshape(dfdu(t,x,u) * d2udQ2_i, [nx*nTq,nTq]); % ((fu_Q -> fQ_u) * u_Q -> fQ_Q) + (f_u * u_QQ -> f_QQ -> fQ_Q) -> fQ_Q
            
            d2fdqdk = d2fdudk(t,x,u); % fk_u
            d2fdqdk = d2fdqdk(fkUseParams,:) * dudq_i; % (fk_u -> fK_u) * u_Q -> fK_Q
            
            val = [val(fkUseParams,opts.UseParams),                  sparse(nx*nTk, nTs), d2fdqdk, sparse(nx*nTk,nTh);
                   sparse(nx*nTs, nTk+nTs+nTq+nTh);
                   spermute132(d2fdqdk, [nx,nTk,nTq], [nx*nTq,nTk]), sparse(nx*nTq, nTs), d2fdq2,  sparse(nx*nTh,nTh);
                   sparse(nx*nTh, nTk+nTs+nTq+nTh)];
               % fT_T
        end
        
        % Modifies d2fdkdx to relate only to the parameters of interest
        function val = d2fdTdxSub(t, x, u)
            val = d2fdkdx(t,x,u); % fx_k
            d2fdqdx = d2fdudx(t,x,u) * dudq(t); % fx_u * u_q -> fx_q
            val = [val(:,opts.UseParams), sparse(nx*nx, nTs), d2fdqdx(:,opts.UseInputControls), sparse(nx*nx,nTh)]; % fx_T
        end
        
        % Modifies d2fdxdk to relate only to the parameters of interest
        function val = d2fdxdTSub(t, x, u)
            val = spermute132(d2fdTdx(t,x,u), [nx,nx,nT], [nx*nT,nx]); % fx_T -> fT_x
        end
    end
end
