function ic = steadystateCurv(m, con, opts)
% Constants
nx = m.nx;
nu = m.nu;
nk = m.nk;
nq = con.nq;
nTk = sum(opts.UseParams);
nTs = sum(opts.UseSeeds);
nTq = sum(opts.UseInputControls);
nTh = sum(opts.UseDoseControls);
nT  = nTk + nTs + nTq + nTh;

uqUseInputControls = linearslicer([nu,nq], true(nu,1), opts.UseInputControls);

% Construct system
[der, jac, eve] = constructSystem();

% Initial conditions
order = 2;
ic = extractICs(m,con,opts,order);

% Integrate [f; dfdT] over time
sol = accumulateOdeFwdSimp(der, jac, 0, inf, ic, con.private.BasalDiscontinuities, 0, 1:nx, opts.RelTol, opts.AbsTol(1:nx+nx*nT+nx*nT*nT), [], eve, @(cum_sol)true);

% Return steady-state value
ic = sol.ye;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f and dfdT and d2fdT2 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, eve] = constructSystem()
        
        IT = speye(nT);
        IT2 = speye(nT*nT);
        TPermuteInd = vec(permute(reshape(1:nx*nT*nT, nx,nT,nT), [1,3,2])); % indexes that permute the last two T of fTT_
        fkUseParams = reshape(1:nx*nk, nx,nk);
        fkUseParams = vec(fkUseParams(:,opts.UseParams)); % indexes that remove inactive k from fk_
%         fqUseControls = reshape(1:nx*nTq, nx,nTq);
%         fqUseControls = vec(fqUseControls(:,opts.UseControls{1})); % indexes that remove inactive q from fq_

        dxdTStart = nx+1;
        dxdTEnd   = nx+nx*nT;
        d2xdT2Start = nx+nx*nT+1;
        d2xdT2End   = nx+nx*nT+nx*nT*nT;
        
%         f       = m.f;
%         dfdx    = m.dfdx;
%         dfdu    = m.dfdu;
%         dfdk    = m.dfdk;
%         dfdT    = @dfdTSub;
%         d2fdx2  = m.d2fdx2;
%         d2fdkdx = m.d2fdkdx;
%         d2fdT2  = @d2fdT2Sub;
%         d2fdudx = m.d2fdudx;
%         d2fdTdx = @d2fdTdxSub;
%         d2fdxdT = @d2fdxdTSub;

        
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
        uf      = con.private.basal_u;
        dudq    = con.private.basal_dudq;
        d2udq2  = con.private.basal_d2udq2;
        
        der = @derivative;
        jac = @jacobian;
        eve = @events;
        
        % Derivative of [x; dxdT; d2xdT2] with respect to time
        function val = derivative(t, joint)
            u = uf(t);            
            x = joint(1:nx);
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % xT_ -> x_T
            d2xdT2 = reshape(joint(d2xdT2Start:d2xdT2End), nx,nT*nT); % xTT_ -> x_TT
            
            % Derivative of dxdT
            % dxdTdot = dfdx *{x.x} dxdT + dxdT
            dxdTdot = dfdx(-1,x,u) * dxdT + dfdT(t,x,u); % f_x * x_T + f_T -> f_T
            
            % Derivative of d2xdT2
            % d2xdT2dot = dfdx *{x.x} d2xdT1dT2 + (d2fdx2 *{x.x} dxdT2) *{x.x} dxdT1 + d2fdT2dx *{x.x} dxdT1 + d2fdT1dx *{x.x} dxdT2 + d2fdT1dT2
            temp = sparse(d2fdxdT(t,x,u) * dxdT); % fT_x * x_T -> fT_T
            temp = temp + spermute132(temp, [nx,nT,nT], [nx*nT,nT]); % fT_T + (fT_T -> fT_T) -> fT_T
            d2xdT2dot = sparse(d2fdx2(-1,x,u) * dxdT); % fx_x * x_T -> fx_T
            d2xdT2dot = spermute132(d2xdT2dot, [nx,nx,nT], [nx*nT,nx]); % fx_T -> fT_x
            d2xdT2dot = sparse(d2xdT2dot * dxdT); %fT_x * x_T -> fT_T
            d2xdT2dot = reshape(dfdx(-1,x,u) * d2xdT2, nx*nT,nT) + d2xdT2dot + temp + d2fdT2(t,x,u); % (f_x * x_TT -> f_TT -> fT_T) + fT_T + fT_T + fT_T -> fT_T
            
            % Combine
            val = [f(-1,x,u); vec(dxdTdot); vec(d2xdT2dot)];
        end
        
        % Jacobian of [x; dxdT; d2xdT2] derivative
        function val = jacobian(t, joint)
            u = uf(t);            
            x = joint(1:nx); % x_
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % x_T
            d2xdT2 = reshape(joint(d2xdT2Start:d2xdT2End), nx,nT*nT); % xTT_ -> x_TT
            
            % Compute d/dx(dfdT)
            d2fdxdTVal = sparse(d2fdx2(-1,x,u) * dxdT) + d2fdTdx(t,x,u); % fx_T
            d2fdxdTVal = spermute132(d2fdxdTVal, [nx,nx,nT], [nx*nT,nx]);
            
            % Compute d/dx(d2fdT2)
            d3fdxdT2 = sparse(nx*nT*nT,nx);% TEMP
            
            % Compute d/dxdT(d2fdT2)
            % d/dxdT(d2fdT2) = (d2fdx2 *{x.x} dxdT2) * IT1 + (d2fdx2 *{x.x} dxdT1) * IT2 + d2fdT2dx * IT1 + d2fdT1dx * IT2
            d3fddxdTT2 = sparse(d2fdx2(-1,x,u) * dxdT); % fx_x * x_T -> fx_T
            d3fddxdTT2 = spermute132(d3fddxdTT2, [nx,nx,nT], [nx*nT,nx]); % fx_T -> fT_x
            d3fddxdTT2 = d3fddxdTT2 + d2fdxdT(t,x,u); % fT_x + fT_x -> fT_x
            d3fddxdTT2 = kron(IT, d3fddxdTT2); % T_T kron fT_x -> fTT_xT
            d3fddxdTT2 = d3fddxdTT2 + d3fddxdTT2(TPermuteInd,:); % fTT_xT + (fTT_xT -> fTT_xT) -> fTT_xT
            
            % Combine
            val = [dfdx(-1,x,u), sparse(nx,nx*nT),      sparse(nx,nx*nT*nT);
                   d2fdxdTVal,  kron(IT, dfdx(-1,x,u)), sparse(nx*nT,nx*nT*nT);
                   d3fdxdT2,    d3fddxdTT2,            kron(IT2, dfdx(-1,x,u))];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = dfdk(-1,x,u); % f_k
            dudq_i = dudq(t); % u_q
            dfdTq = dfdu(-1,x,u) * dudq_i(:,opts.UseInputControls); % f_u * (u_q -> u_Q) -> f_Q
            val = [val(:,opts.UseParams), sparse(nx,nTs), dfdTq, sparse(nx,nTh)]; % f_T
        end
        
        % Modifies d2fdkdk to relate only to the parameters of interest
        function val = d2fdT2Sub(t, x, u)
            val = d2fdk2(-1,x,u); % fk_k
            
            dudq_i = dudq(t); % u_q
            dudq_i = dudq_i(:,opts.UseInputControls); % u_q -> u_Q
            
            d2udQ2_i = d2udq2(t); % uq_q
            d2udQ2_i = reshape(d2udQ2_i(uqUseInputControls,opts.UseInputControls), [nu,nTq*nTq]); % uq_q -> uQ_Q -> u_QQ
            
            d2fdq2 = d2fdu2(-1,x,u) * dudq_i; % fu_u * u_Q -> fu_Q
            d2fdq2 = spermute132(d2fdq2, [nx,nu,nTq], [nx*nTq,nu]) * dudq_i + reshape(dfdu(-1,x,u) * d2udQ2_i, [nx*nTq,nTq]); % ((fu_Q -> fQ_u) * u_Q -> fQ_Q) + (f_u * u_QQ -> f_QQ -> fQ_Q) -> fQ_Q
            
            d2fdqdk = d2fdudk(-1,x,u); % fk_u
            d2fdqdk = d2fdqdk(fkUseParams,:) * dudq_i; % (fk_u -> fK_u) * u_Q -> fK_Q
            
            val = [val(fkUseParams,opts.UseParams),                  sparse(nx*nTk, nTs), d2fdqdk, sparse(nx*nTk,nTh);
                   sparse(nx*nTs, nTk+nTs+nTq+nTh);
                   spermute132(d2fdqdk, [nx,nTk,nTq], [nx*nTq,nTk]), sparse(nx*nTq, nTs), d2fdq2,  sparse(nx*nTh,nTh);
                   sparse(nx*nTh, nTk+nTs+nTq+nTh)];
               % fT_T
        end
        
        % Modifies d2fdkdx to relate only to the parameters of interest
        function val = d2fdTdxSub(t, x, u)
            val = d2fdkdx(-1,x,u); % fx_k
            d2fdqdx = d2fdudx(-1,x,u) * dudq(t); % fx_u * u_q -> fx_q
            val = [val(:,opts.UseParams), sparse(nx*nx, nTs), d2fdqdx(:,opts.UseInputControls), sparse(nx*nx,nTh)]; % fx_T
        end
        
        % Modifies d2fdxdk to relate only to the parameters of interest
        function val = d2fdxdTSub(t, x, u)
            val = spermute132(d2fdTdx(t,x,u), [nx,nx,nT], [nx*nT,nx]); % fx_T -> fT_x
        end
        
        % Steady-state event
        function [value, isTerminal, direction] = events(t, joint)
            u = uf(t);
            x = joint(1:nx); % x_

            % Absolute change
            absDiff = con.private.TimeScale * f(-1,x,u); % Change over an entire simulation
            
            % Relative change
            relDiff = absDiff ./ x;
            
            % Either absolute change or relative change must be less than
            % the tolerance for all species
            value = max(min(abs(absDiff) - opts.AbsTol(1:nx), abs(relDiff) - opts.RelTol));
            if value < 0
                value = 0;
            end
            
            % Always end and only care about drops below the threshold
            isTerminal = true;
            direction = -1;
        end
    end
end
