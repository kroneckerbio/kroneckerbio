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

if opts.Verbose; fprintf('Integrating sensitivities:\n'); end
for iCon = 1:nCon
    if verboseAll; tic; end
    
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

    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.ObjWeights = opts.ObjWeights(:,iCon);
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
    end
    discD = zeros(inT,nObj*nDisc); % T_ot
    for iObj = 1:nObj
        for iDisc = 1:nDisc
            objDiscD = vec(vec(obj(iObj,iCon).dGdx(discreteTimes(iDisc), sol)).' * reshape(dxdT(:,iDisc), nx, inT)); % _x * x_T --> _T --> T_
            temp = vec(obj(iObj,iCon).dGdk(discreteTimes(iDisc), sol)); % k_ % partial dGdk(i)
            temp = [temp(opts.UseParams); zeros(inTx+inTq,1)]; % k_ --> T_
            objDiscD = objDiscD + temp; % T_
            discD(:,(iObj-1)*nDisc + iDisc) = opts.ObjWeights(iObj,iCon) * objDiscD; % T_ as a row in T_ot
        end
    end
    
    % Sorting by abs value to minimize numerical error
    [unused I] = sort(abs(discD),2);
    for iT = 1:inT
        discD(iT,:) = discD(iT,I(iT,:));
    end
    discD = sum(discD,2);
    
    % Sum discrete and continuous terms
    % Rate parameters are the same
    D(1:nTk) = D(1:nTk) + contD(1:nTk) + discD(1:nTk);
    
    if opts.UseModelICs
        D(nTk+1:nTk+nTx) = D(nTk+1:nTk+nTx) + contD(nTk+1:nTk+inTx) + discD(nTk+1:nTk+inTx);
    else
        % x0 parameters are different
        D(Txind+1:Txind+inTx) = contD(nTk+1:nTk+inTx) + discD(nTk+1:nTk+inTx);
        % Increment index of variable ICs
        Txind = Txind + inTx;
    end
    
    if opts.UseModelInputs
        D(nTk+nTx+1:end) = D(nTk+nTx+1:end) + contD(nTk+inTx+1:end) + discD(nTk+inTx+1:end);
    else
        % q parameters are different
        D(Tqind+1:Tqind+inTq) = contD(nTk+inTx+1:end) + discD(nTk+inTx+1:end);
        Tqind = Tqind + inTq;
    end
    
    % Return the adaptive abstol, if requested
%     if nargout >= 3
%         abstol{iCon} = opts.RelTol * abstolObjSimple(m, con(iCon), obj(:,iCon), sol);
%     end
    
    if verboseAll; fprintf('iCon = %d\t|dGdT| = %g\tTime = %0.2f\n', iCon, norm(contD + discD), toc); end
end

if opts.Verbose; fprintf('Summary: |dGdT| = %g\n', norm(D)); end