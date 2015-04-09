function [G, D] = computeObjSens(m, con, obj, opts)
% Switch to Adjoint method if requested
if opts.UseAdjoint
    [G, D] = computeObjSensAdj(m, con, obj, opts);
    return
end

% Continue with forward method
verboseAll = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nTk = nnz(opts.UseParams);
nTs = nnz(opts.UseSeeds);
nTq = nnz(cat(1,opts.UseInputControls{:}));
nTh = nnz(cat(1,opts.UseDoseControls{:}));
nT  = nTk + nTs + nTq + nTh;
nCon = numel(con);
nObj = size(obj,1);

% Initialize variables
G = 0;
D = zeros(nT,1);

if opts.Verbose; fprintf('Integrating sensitivities:\n'); end
for iCon = 1:nCon
    if verboseAll; tic; end
    
    % Modify opts structure
    intOpts = opts;
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.ObjWeights = opts.ObjWeights(:,iCon);
    
    UseSeeds_i = opts.UseSeeds(:,iCon);
    intOpts.UseSeeds = UseSeeds_i;
    inTs = nnz(UseSeeds_i);
    
    UseInputControls_i = opts.UseInputControls{iCon};
    intOpts.UseInputControls = UseInputControls_i;
    inTq = nnz(UseInputControls_i);
    
    UseDoseControls_i = opts.UseDoseControls{iCon};
    intOpts.UseDoseControls = UseDoseControls_i;
    inTh = nnz(UseDoseControls_i);
    
    inT = nTk + inTs + inTq + inTh;

    tGet = opts.tGet{iCon};
    
    % Integrate
    if opts.continuous(iCon) && opts.complex(iCon)
        sol = integrateObjSens(m, con(iCon), obj(:,iCon), intOpts);
    elseif opts.complex(iCon)
        sol = integrateSens(m, con(iCon), intOpts);
    elseif opts.continuous(iCon)
        sol = integrateObjSensSelect(m, con(iCon), obj(:,iCon), tGet, intOpts);
    else
        sol = integrateSensSelect(m, con(iCon), tGet, intOpts);
    end
    
    % Add fields for prior objectives
    sol.UseParams = opts.UseParams;
    sol.UseSeeds = UseSeeds_i;
    sol.UseInputControls = UseInputControls_i;
    sol.UseDoseControls = UseDoseControls_i;
    
    % *Compute G*
    % Extract continuous term
    if opts.continuous(iCon)
        contG = sol.y(nx+1,end);
    else
        contG = 0;
    end
    
    % Compute discrete term
    discG = 0;
    discreteTimes = zeros(0,1);
    for iObj = 1:nObj
        [iDiscG, temp] = obj(iObj,iCon).G(sol);
        discreteTimes = [discreteTimes; vec(temp)];
        discG = discG + opts.ObjWeights(iObj,iCon) * iDiscG;
    end

    % Remove repetitive discreteTimes
    discreteTimes = unique(discreteTimes);
    nDisc = numel(discreteTimes);

    % Add to cumulative goal value
    G = G + contG + discG;
    
    % *Compute D*
    % Extract continuous term
    if opts.continuous(iCon)
        dGdTStart = nx+1+nx*inT+1;
        dGdTEnd   = nx+1+nx*inT+inT;
        contD = sol.y(dGdTStart:dGdTEnd,end);
    else
        contD = zeros(nT,1);
    end
    
    % Compute discrete term
    dxdTStart = nx+opts.continuous(iCon)+1;
    dxdTEnd   = nx+opts.continuous(iCon)+nx*inT;
    if ~isempty(discreteTimes) % deval crashes on empty t
        dxdT = deval(sol, discreteTimes, dxdTStart:dxdTEnd); % xT_t
    else
        dxdT = zeros(nx*nT,0);
    end
    discD = zeros(inT,nObj*nDisc); % T_ot
    for iObj = 1:nObj
        for iDisc = 1:nDisc
            objDiscD = vec(vec(obj(iObj,iCon).dGdx(discreteTimes(iDisc), sol)).' * reshape(dxdT(:,iDisc), nx, inT)); % _x * x_T --> _T --> T_
            dGdk = obj(iObj,iCon).dGdk(discreteTimes(iDisc), sol); % k_ % partial dGdk(i)
            dGds = obj(iObj,iCon).dGds(discreteTimes(iDisc), sol); % s_ % partial dGds(i)
            dGdq = obj(iObj,iCon).dGdq(discreteTimes(iDisc), sol); % q_ % partial dGdq(i)
            dGdh = obj(iObj,iCon).dGdh(discreteTimes(iDisc), sol); % h_ % partial dGdh(i)
            objDiscD = objDiscD + [dGdk(opts.UseParams); dGds(UseSeeds_i); dGdq(UseInputControls_i); dGdh(UseDoseControls_i)]; % T_
            discD(:,(iObj-1)*nDisc + iDisc) = opts.ObjWeights(iObj,iCon) * objDiscD; % T_ as a row in T_ot
        end
    end
    
    % Sorting by abs value to minimize numerical error
    [unused, I] = sort(abs(discD),2);
    for iT = 1:inT
        discD(iT,:) = discD(iT,I(iT,:));
    end
    discD = sum(discD,2);
    
    % Sum discrete and continuous terms
    Di = zeros(nT,1);
    
    % Rate parameters are the same
    Di(1:nTk) = Di(1:nTk) + contD(1:nTk) + discD(1:nTk);
    
    % s parameters are different
    Tsind = nTk + sum(sum(opts.UseSeeds(:,1:iCon-1)));
    Di(Tsind+1:Tsind+inTs) = contD(nTk+1:nTk+inTs) + discD(nTk+1:nTk+inTs);
    
    % q parameters are different
    Tqind = nTk + nTs + sum(cellfun(@sum, opts.UseInputControls(1:iCon-1)));
    Di(Tqind+1:Tqind+inTq) = contD(nTk+inTs+1:nTk+inTs+inTq) + discD(nTk+inTs+1:nTk+inTs+inTq);
    
    % h parameters are different
    Thind = nTk + nTs + nTq + sum(cellfun(@sum, opts.UseDoseControls(1:iCon-1)));
    Di(Thind+1:Thind+inTh) = contD(nTk+inTs+inTq+1:end) + discD(nTk+inTs+inTq+1:end);
    
    % Update sum
    D = D + Di;
    
    if verboseAll; fprintf('iCon = %d\t|dGdT| = %g\tTime = %0.2f\n', iCon, norm(contD + discD), toc); end
end

if opts.Verbose; fprintf('Summary: |dGdT| = %g\n', norm(D)); end
