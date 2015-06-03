function ic = extractICs(m,con,opts,order)
% function ic = extractICs(m,con,opts,order) extracts the initial
% conditions and their derivatives with respect to fitted rate, seed,
% input, and dosing parameters, according to the specified derivative
% order. opts need only contain the UseParams, UseSeeds, UseInputControls,
% and UseDoseControls fields to indicate which parameters are to be fitted.
% If order == 0, these fields are not required.

% Initial conditions
x0 = m.x0(con.s);

if order >= 1
    
    % Constants
    nx = m.nx;
    ns = m.ns;
    nTk = sum(opts.UseParams);
    nTs = sum(opts.UseSeeds);
    nTq = sum(opts.UseInputControls);
    nTh = sum(opts.UseDoseControls);
    nT  = nTk + nTs + nTq + nTh;
    
    % Initial effect of rates on sensitivities is 0
    dx0dTk = zeros(nx, nTk); % Active rate parameters
    
    % Initial effect of seeds on states is dx0ds
    dx0dTs = m.dx0ds(con.s);
    dx0dTs = dx0dTs(:,opts.UseSeeds);
    
    % Initial effect of qs on sensitivities is 0
    dx0dTq = zeros(nx, nTq);
    
    % Initial effect of hs on sensitivities is 0
    dx0dTh = zeros(nx, nTh);
    
    if order >= 2
        
        % Initial curvatures are zero, except for seeds
        
        % Get indices corresponding to fitted seeds in T
        T_seedindices = nTk+1:nTk+nTs;
        
        % Determine which elements of d2x0dT2 involve only seeds
        d2x0dT2_isseeds = false(nx,nT,nT);
        d2x0dT2_isseeds(1:nx,T_seedindices,T_seedindices) = true;
        
        % Get derivatives of x0 with respect to all seeds
        d2x0ds2 = m.d2x0ds2(con.s);
        
        % Filter out non-fitted seeds
        d2x0dT2_s = reshape(d2x0ds2,nx,ns*ns);
        UseSeeds_ss = false(ns,ns);
        UseSeeds_ss(opts.UseSeeds,opts.UseSeeds) = true;
        d2x0dT2_s = d2x0dT2_s(:,UseSeeds_ss(:));
        d2x0dT2_s = reshape(d2x0dT2_s,nx*nTs,nTs);
        
        % Set nonzero elements of d2x0dT2
        d2x0dT2 = zeros(nx*nT*nT,1);
        d2x0dT2(vec(d2x0dT2_isseeds)) = vec(d2x0dT2_s);
        
    else
        
        d2x0dT2 = [];
        
    end
    
else
    
    dx0dTk = [];
    dx0dTs = [];
    dx0dTq = [];
    dx0dTh = [];
    d2x0dT2 = [];
    
end

% Combine them into a vector
ic = [x0; vec([dx0dTk, dx0dTs, dx0dTq, dx0dTh]); d2x0dT2];