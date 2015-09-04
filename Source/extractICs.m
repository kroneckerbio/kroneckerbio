function ic = extractICs(m,con,opts,order)
% function ic = extractICs(m,con,opts,order) extracts the initial
% conditions and their derivatives with respect to fitted rate, seed,
% input, and dosing parameters, according to the specified derivative
% order. opts need only contain the UseParams, UseSeeds, UseInputControls,
% and UseDoseControls fields to indicate which parameters are to be fitted.
% If order == 0, these fields are not required.

s = con.s;

% Initial conditions
x0 = m.x0(s);

if order >= 1
    % Constants
    nx = m.nx;
    ns = m.ns;
    nk = m.nk;
    
    UseParams = opts.UseParams;
    UseSeeds = opts.UseSeeds;
    
    nTk = sum(UseParams);
    nTs = sum(UseSeeds);
    nTq = sum(opts.UseInputControls);
    nTh = sum(opts.UseDoseControls);
    nT  = nTk + nTs + nTq + nTh;
    
    % Initial effect of rates on states is dx0dk
    dx0dk = m.dx0dk(s);
    dx0dTk = dx0dk(:,UseParams);
    
    % Initial effect of seeds on states is dx0ds
    dx0ds = m.dx0ds(s);
    dx0dTs = dx0ds(:,UseSeeds);
    
    % Initial effect of qs on states is 0
    dx0dTq = zeros(nx, nTq);
    
    % Initial effect of hs on states is 0
    dx0dTh = zeros(nx, nTh);
    
    dx0dT = [dx0dTk, dx0dTs, dx0dTq, dx0dTh];
    
    if order >= 2
        d2x0ds2 = reshape(full(m.d2x0ds2(s)), nx,ns,ns);
        d2x0dTs2 = d2x0ds2(:,UseSeeds,UseSeeds);
        
        d2x0dk2 = reshape(full(m.d2x0dk2(s)), nx,nk,nk);
        d2x0dTk2 = d2x0dk2(:,UseParams,UseParams);
        
        d2x0dkds = reshape(full(m.d2x0dkds(s)), nx,ns,nk);
        d2x0dTkdTs = d2x0dkds(:,UseSeeds,UseParams);
        
        d2x0dsdk = reshape(full(m.d2x0dsdk(s)), nx,nk,ns);
        d2x0dTsdTk = d2x0dsdk(:,UseParams,UseSeeds);
        
        d2x0dT2 = zeros(nx,nT,nT);
        d2x0dT2(:,1:nTk,1:nTk) = d2x0dTk2;
        d2x0dT2(:,nTk+(1:nTs),nTk+(1:nTs)) = d2x0dTs2;
        d2x0dT2(:,nTk+(1:nTs),1:nTk) = d2x0dTkdTs;
        d2x0dT2(:,1:nTk,nTk+(1:nTs)) = d2x0dTsdTk;
    else
        d2x0dT2 = [];
    end
else
    dx0dT = [];
    d2x0dT2 = [];
end

% Combine them into a vector
ic = [x0; vec(dx0dT); vec(d2x0dT2)];
