function int = integrateLnaComp(m, con, tF, eve, fin, opts)

% Constants
nr = m.nr;
nx = m.nx;
ny = m.ny;

y = m.y;
u = con.u;
dydx = m.dydx;

upperVInd = upperInd(nx);
lowerVInd = lowerInd(nx);
nVx = numel(upperVInd);

% Construct system
[der, jac, del] = constructSystem();

% V is symmetric. Only integrate the upper half of the matrix
if ~con.SteadyState
    order = 0;
    ic = extractICs(m,con,opts,order);
    
    V0 = opts.V0(upperVInd);

    ic = [ic; V0];
else
    error('Steady-state not yet implemented for LNA')
    ic = steadystateSys(m, con, opts);
end

% Integrate x over time
sol = accumulateOdeFwdComp(der, jac, 0, tF, ic, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+numel(upperVInd)), del, eve, fin);

% Work down
int.Type = 'Integration.Lna.Complex';
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

int.t = sol.x;
int.x = @(t)devals(sol, t);
int.u = con.u;
int.y = @(t)y(t, devals(sol, t), u(t));

int.Vx = @(t)evaluate_state_variance(sol,t);
int.Vy = @(t)evaluate_output_variance(sol,t);

int.ie = sol.ie;
int.te = sol.xe;
int.xe = sol.ye(1:nx,:);
int.ue = u(int.te);
int.ye = y(int.te, int.xe, int.ue);

int.sol = sol;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating the linear noise approximation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        S       = m.S;
        r       = m.r;
        drdx    = m.drdx;
        f       = m.f;
        dfdx    = m.dfdx;
        d2fdx2  = m.d2fdx2;
        uf      = con.u;
        d       = con.d;
        x0      = m.x0;
        nd      = m.ns;
        
        Ix = speye(nx);
        
        VStart = nx+1;
        VEnd   = nx+numel(upperVInd);
        diagVInd = sub2ind([nx,nx], vec(1:nx), vec(1:nx));
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [x;V] with respect to time
        function val = derivative(t, joint)
            x = joint(1:nx); % x_
            V = zeros(nx,nx);
            V(upperVInd) = joint(VStart:VEnd); % x_x
            V(lowerVInd) = joint(VStart:VEnd); % x_x
            u_t = uf(t); % u_
            
            % xdot = f
            % Vdot = dfdx *{x.x} V +{f+V.x;V.x+f} dfdx *{x.x} V +{f+x;f+x} S *{rxr} r *{r.r} S
            Vdot = dfdx(t,x,u_t) * V; % f_x * v1_V.x -> f_V.x
            Vdot = Vdot + Vdot.' + bsxfun(@times, S, r(t,x,u_t).') * S.'; % f_V.x + V.x_f + x_r .* r * r_x -> x_x
            
            val = [f(t,x,u_t); Vdot(upperVInd)];
        end
        
        % Jacobian of x derivative
        function val = jacobian(t, joint)
            x = joint(1:nx); % x_
            V = zeros(nx,nx);
            V(upperVInd) = joint(VStart:VEnd); % x_x
            V(lowerVInd) = joint(VStart:VEnd); % x_x
            u_t = uf(t); % u_
            
            % dxdot/dx = dfdx
            % dxdot/dV = 0
            % dVdot/dx = d2fdx2 *{x1.x} V +{f+V.x;V.x+f;x2+x2} d2fdx2 *{x1.x} V +{f+S.x;f+S.x;x2+drdx.x} S *{rxr} drdx *{r.r} S
            % dVdot/dV = dfdx *{x.x} Iv +{f+V1.x;V1.x+f} dfdx *{x.x} dfdx
            dVdotdx = d2fdx2(t,x,u_t) * V; % fx_x * V.x_V.x -> fx_V.x
            dVdotdx = spermute132(dVdotdx, [nx,nx,nx], [nx*nx,nx]) + reshape(dVdotdx.', [nx*nx,nx]); % (fx_V.x -> f,V.x_x) + (fx_V.x -> V.x_fx -> V.x,f_x) -> ff_x
            dVdotdx = dVdotdx + reshape(S * reshape(bsxfun(@times, vec(S.'), repmat(drdx(t,x,u_t), nx,1)), nr,nx*nx), nx*nx,nx); % ff_x + (f_r * ((f_r -> r_f -> rf_) .* (r_x -> rf_x) -> rf_x -> r_fx) -> f_fx -> ff_x) -> ff_x
            dVdotdx = dVdotdx(upperVInd,:);
            
            dVdordV = kron(Ix, dfdx(t,x,u_t)) + kron(dfdx(t,x,u_t), Ix) + sparse(repmat(vec(1:nx*nx), nx,1), vec(repmat(1:nx*nx, nx,1)), vec(repmat(dfdx(t,x,u_t), nx,1)), nx*nx,nx*nx) + sparse(vec(repmat(1:nx*nx, nx,1)), repmat(vec(1:nx*nx), nx,1), vec(repmat(dfdx(t,x,u_t).', nx,1)), nx*nx,nx*nx);
            dVdordV(:,diagVInd) = dVdordV(:,diagVInd) ./ 2; % Diagonal does not get copied, so it has half the normal effect
            dVdordV = dVdordV(upperVInd,upperVInd);
            
            val = [dfdx(t,x,u_t), zeros(nx, numel(upperVInd));
                   dVdotdx,     dVdordV];
        end
        
        % Dosing
        function val = delta(t, joint)
            val = [x0(d(t)) - x0(zeros(nd,1)); zeros(numel(upperVInd),1)];
        end
        
    end

    function val = evaluate_state_variance(sol, t)
        nt = numel(t);
        
        extract = devals(sol, t,nx+1:nx+nVx);
        val = zeros(nx*nx,nt);
        val(upperVInd,:) = extract;
        val(lowerVInd,:) = extract; % xx_t
    end

    function val = evaluate_output_variance(sol, t)
        nt = numel(t);
        
        xs = devals(sol, t,1:nx);
        Vxs = evaluate_state_variance(sol, t);
        
        val = zeros(ny*ny,nt);
        for i = 1:nt
            x_i = xs(:,i); % x_
            u_i =  u(t(i)); % u_
            Vx_i = reshape(Vxs(:,i), nx,nx); % xx_ -> x_x
            dydx_i = dydx(t(i), x_i, u_i);
            val(:,i) = vec(dydx_i * Vx_i * dydx_i.'); % y_x * x_x * (y_x -> x_y) -> y_y -> yy_
        end
    end
end
