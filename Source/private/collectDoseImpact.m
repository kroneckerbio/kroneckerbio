function ic = collectDoseImpact(m, con, t, order, UseParams, UseSeeds, UseInputControls, UseDoseControls)
% Collect the impact on the state of the system given an applied dose

if nargin < 8
    UseDoseControls = [];
    if nargin < 7
        UseInputControls = [];
        if nargin < 6
            UseSeeds = [];
            if nargin < 5
                UseParams = [];
            end
        end
    end
end

% Constants
nx = m.nx;
ns = m.ns;
nk = m.nk;
nh = con.nh;

% Dose effect
d = con.d(t);
d0 = zeros(ns,1);
x0 = m.x0(d) - m.x0(d0);

if order >= 1
    nTk = sum(UseParams);
    nTs = sum(UseSeeds);
    nTq = sum(UseInputControls);
    nTh = sum(UseDoseControls);
    nT  = nTk + nTs + nTq + nTh;
    
    % Initial effect of rates on states is dx0dk
    dx0dk = m.dx0dk(d) - m.dx0dk(d0);
    dx0dTk = dx0dk(:,UseParams);
    
    % Initial effect of seeds on states is 0
    dx0dTs = zeros(nx, nTs);
    
    % Initial effect of qs on states is 0
    dx0dTq = zeros(nx, nTq);
    
    % Initial effect of hs on states is dx0ds * dddh
    dx0ds = m.dx0ds(d);
    dddh = con.dddh(t);
    dx0dTh = dx0ds * dddh(:,UseDoseControls);
    
    dx0dT = [dx0dTk, dx0dTs, dx0dTq, dx0dTh];
    
    if order >= 2
        inds_k = 1:nTk;
        inds_h = nTk+nTs+nTq+(1:nTh);
        dhUseDoseControls = linearslicer([ns,nh], true(ns,1), UseDoseControls);

        d2x0dk2 = m.d2x0dk2(d) - m.d2x0dk2(d0);
        d2x0dTk2 = reshape(full(d2x0dk2), nx,nk,nk);
        d2x0dTk2 = d2x0dTk2(:,UseParams,UseParams);
        
        % d2xdh2 = (dxdd2dd1 *{d.d} dd2dh2) *{d.d} dddh1 + dxdd *{d.d} d2ddh2dh1
        d2x0ds2 = m.d2x0ds2(d); % xs_s
        dddh2 = con.d2ddh2(t); % sh_h
        d2x0dTh2 = d2x0ds2 * dddh; % xs_s * s_h -> xs_h
        d2x0dTh2 = spermute132(d2x0dTh2(:,UseDoseControls), [nx,ns,nTh], [nx*nTh,ns]); % xs_h -> xs_H -> xH_s
        d2x0dTh2 = d2x0dTh2 * dddh; % xH_s * s_h -> xH_h
        d2x0dTh2 = reshape(full(d2x0dTh2(:,UseDoseControls)), nx,nTh,nTh) ... % xH_h -> xH_H -> x_H_H
            + reshape(full(dx0ds * reshape(dddh2(dhUseDoseControls, UseDoseControls), ns,nTh*nTh)), nx,nTh,nTh); % x_s * (sh_h -> sH_H -> s_HH) -> x_HH -> x_H_H
        
        d2x0dsdk = m.d2x0dsdk(d);
        d2x0dThdTk = d2x0dsdk * dddh; % xk_s * s_h -> xk_h
        d2x0dThdTk = reshape(full(d2x0dThdTk), nx,nk,nh); % xk_h -> x_k_h
        d2x0dThdTk = d2x0dThdTk(:,UseParams,UseDoseControls); % x_k_h -> x_K_H
        
        d2x0dT2 = zeros(nx,nT,nT);
        d2x0dT2(:,inds_k,inds_k) = d2x0dTk2;
        d2x0dT2(:,inds_h,inds_h) = d2x0dTh2;
        d2x0dT2(:,inds_h,inds_k) = permute(d2x0dThdTk, [1,3,2]); % x_K_H -> x_H_K
        d2x0dT2(:,inds_k,inds_h) = d2x0dThdTk;
    else
        d2x0dT2 = [];
    end
else
    dx0dT = [];
    d2x0dT2 = [];
end

% Combine them into a vector
ic = [x0; vec(dx0dT); vec(d2x0dT2)];
