function [optimalGRatio optimalDRatio] = integrateOptimalAbsTolSensSimple(m, con, obj, opts)
%% Constants
nx = m.nx;
nTk = sum(opts.UseParams);
nTx = sum(sum(opts.UseICs));
nT  = nTk + nTx;
nObj = numel(obj);

verbose = opts.Verbose;
opts.Verbose = max(opts.Verbose-1,0);

absTolStart = nx+nx*nT+1;
absTolEnd   = nx+nx*nT+nx*nx+nx*nT*nx;
dxdTStart   = nx+1;
dxdTEnd     = nx+nx*nT;
dxdx0Start  = 1;
dxdx0End    = nx*nx;
dsdx0Start  = nx*nx+1;
dsdx0End    = nx*nx+nx*nT*nx;

%% Integrate the basic system
% Complete integration required even if objective function doesn't need it
if verbose; fprintf('Integrating initial sensitivities...'); end
sol = integrateSens(m, con, opts);
if verbose; fprintf('done.\n'); end

nt = numel(sol.x) - 1; % First point has no error

%% Chop down coverage of integration
if opts.Coverage <= 1
    nPick = ceil(opts.Coverage * nt);
elseif opts.Coverage > 1
    nPick = ceil(opts.Coverage);
end
selectIndex = ceil((0:nPick-1) * nt / nPick) + 2; % Add 1 to undo 0 index, add 1 more to skip first point

if verbose; fprintf('Covering %d out of %d times...\n', nPick, nt); end

%% Objective sensitivity to species at all times
% Discrete Times
discreteTimes = [];
for iObj = 1:nObj
    discreteTimes = [discreteTimes; vec(obj(iObj).DiscreteTimes)];
end
discreteTimes = unique(discreteTimes);
nDisc = numel(discreteTimes);

% Precompute dG/dx at all times
dGdx = zeros(nx,nDisc);
for iDisc = 1:nDisc
    for iObj = 1:nObj
        dGdx(:,iDisc) = dGdx(:,iDisc) + opts.ObjWeights(iObj) * vec(obj(iObj).dGdx(discreteTimes(iDisc), sol));
    end
end

%% Sensitivities at all timepoints
effect = zeros(nx+nT*nx,nPick);

% Initial conditions are the same for all runs
ic = zeros(nx*nx+nx*nx*nT,1);
ic(dxdx0Start:dxdx0End) = speye(nx);

% Input
if opts.UseModelInputs
    u = m.u;
else
    u = con.u;
end

for it = 1:nPick
    % Construct system
    [der, jac] = constructSystem();
    
    % Set the integration intervals
    t0 = sol.x(selectIndex(it));
    Disc0 = find(discreteTimes >= t0, 1); % First index in discreteTimes that matters
    
    % Integrate [dx/dx0; d2x/dx0dv]
    if verbose; fprintf('Integrating from time %d...', t0); end
    % There is actual integration to be done
    soli = accumulateOde(der, jac, t0, con.tF, ic, u, con.Discontinuities, [], opts.RelTol, opts.AbsTol(absTolStart:absTolEnd));
    if verbose; fprintf('done.\n'); end
    
    % Compute the effect over remaining time
    for iu = Disc0:nDisc % index into discreteTimes
        effect(:,it) = effect(:,it) + vec(dGdx(:,iu).' * reshape(soli.y(:,iu-Disc0+1), nx,nx+nT*nx)); % multiply by dG/dx, and sum over time
    end
end

% The effect is absolute
effect = abs(effect);

clear soli % very large matrix

%% Maximum impact
maximpact = zeros(1+nT,1);
for it = 1:nPick
    % The impact is absolute
    x = abs(sol.y(1:nx,selectIndex(it))); % x_
    xbox = spdiags(x,0,nx,nx); % x_x diagonal
    dxdT = abs(reshape(sol.y(nx+1:nx+nx*nT,selectIndex(it)),nx,nT)); % x_T
    
    xeffect = effect(1:nx,it);
    xeffectbox = spdiags(xeffect,0,nx,nx);
    dxdpeffect = reshape(effect(nx+1:nx+nx*nT),nT,nx);
    
    % Impact that x has on G
    maximpact(1) = max([maximpact(1); xeffect .* x]); % x_ .* x_ --> x_
    
    % Impact that x has on dG/v
    maximpact(2:end) = max(maximpact(2:end), max(dxdpeffect*xbox, [],2)); % (vx_ --> T_x) * x_x --> T_x --> T_
    
    % Impact that dx/dv has on dG/dv
    maximpact(2:end) = max(maximpact(2:end), vec(max(xeffectbox*dxdT, [],1))); % x_x * x_T --> x_T --> _T --> T_
end

%% Optimal impact from ratio
% Take the maximum over time beforehand
effect = max(effect, [], 2);

% Impact ratio for x's impact on G
optimalGRatio = maximpact(1) ./ effect(1:nx);

% Initialize D
optimalDRatio = zeros(nx+nx*nT,1);

% Impact ratio for x's impact on dG/dv
% There can only be one abstol, take the most conservative
optimalDRatio(1:nx) = min(optimalGRatio, min(spdiags(maximpact(2:end), 0,nT,nT) * reshape(effect(nx+1:nx+nT*nx), nT,nx).^-1, [],1).');

% Impact ratio for dx/dp's impact on dG/dp
% There is no need to take the maximum over v because the entries are
% unique along v. dx/dv only affects dG/dv of the same v
optimalDRatio(nx+1:nx+nx*nT) = (maximpact(2:end) * (effect(1:nx).^(-1)).').'; % T_ * (x_ --> _x) --> T_x --> x_T --> xT_

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [jDer, jJac] = constructSystem()
        
        Ix = speye(nx);
        Ipx = speye(nT*nx);
    
        dfdx    = m.dfdx;
        d2fdx2  = m.d2fdx2;
        d2fdTdx = @d2fdTdxSub;
        
        jDer = @der;
        jJac = @jac;
        
        % Derivative of [dx/dx0; d2x/dx0dv] with respect to time
        function val = der(t, joint, u)
            u = u(t);
            jointold = deval(sol, t, 1:dxdTEnd); % x+xT_
            x = jointold(1:nx); % x_
            dxdT = reshape(jointold(dxdTStart:dxdTEnd), nx,nT); % x_T
            dxdx0 = reshape(joint(dxdx0Start:dxdx0End), nx,nx); % x_x0
            dsdx0 = reshape(joint(dsdx0Start:dsdx0End), nx,nT*nx); % x_Tx0
            
            % Derivative of dx/dx0
            dxdx0dot = vec(dfdx(t,x,u) * dxdx0); % x_x * x_x0 --> x_x0 --> xx0_
            
            % Derivative of ds/dx0
            dsdx0dot = d2fdx2(t,x,u) * dxdT + d2fdTdx(t,x,u); % fx_x * x_T + fx_T --> fx_T
            dsdx0dot = spermute132(dsdx0dot, [nx,nx,nT], [nx*nT,nx]); % fx_T --> fT_x
            dsdx0dot = vec(dsdx0dot * dxdx0) + vec(dfdx(t,x,u) * dsdx0); % (fT_x * x_x0 --> fT_x0 --> fTx0_) + (f_x * x_Tx0 --> f_Tx0 --> fTx0_) --> fTx0_
            
            % Combine all three
            val = [dxdx0dot; dsdx0dot];
        end
        
        % Jacobian of [dx/dx0; d2x/dx0dv] derivative
        function val = jac(t, joint, u)
            u = u(t);
            jointold = deval(sol, t, 1:nx+nx*nT); % x+xT_
            x = jointold(1:nx); % x_
            dxdT = reshape(jointold(dxdTStart:dxdTEnd), nx,nT); % x_T
            
            % Jacobian of ds/dx0 wrt dx/dx0
            d2sdxdx0dx0 = d2fdx2(t,x,u) * dxdT + d2fdTdx(t,x,u); % fx_x * x_T + fx_T --> fx_T
            d2sdxdx0dx0 = spermute132(d2sdxdx0dx0, [nx,nx,nT], [nx*nT,nx]); % fx_T --> fT_x
            
            % Combine all
            val = [kron(Ix, dfdx(t,x,u)), sparse(nx*nx,nx*nT*nx);
                   kron(Ix, d2sdxdx0dx0), kron(Ipx, dfdx(t,x,u))];
        end
        
        % Modifies d2fdkdx to relate only to the parameters of interest
        function val = d2fdTdxSub(t, x, u)
            val = m.d2fdkdx(t,x,u);
            val = [val(:, opts.UseParams) zeros(nx*nx, nTx)];
        end
        
    end
end
