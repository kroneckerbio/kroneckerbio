function G = computeObj(m, con, obj, opts)
% G = computeObj(m, con, obj, opts)
% This function computes the total objective function value for a vector of
% con and a matrix of obj.
verbose = logical(opts.Verbose);
verboseAll = max(verbose-1,0);

% Constants
nx = m.nx;
nTk = sum(opts.UseParams);
nTs = sum(sum(opts.UseSeeds));
nTq = sum(cat(1,opts.UseControls{:}));
nT  = nTk + nTs + nTq;
nCon = numel(con);
nObj = size(obj,1);

% Initialize variables
G = 0;

if verbose; disp('Integrating forward...'); end
for iCon = 1:nCon
    if verboseAll; tic; end
    
    intOpts = opts;
    
    % If opts.UseModelSeeds is false, the number of variables can change
    if opts.UseModelSeeds
        UseSeeds_i = opts.UseSeeds;
    else
        UseSeeds_i = opts.UseSeeds(:,iCon);
    end
    intOpts.UseSeeds = UseSeeds_i;
    inTs = nnz(UseSeeds_i);
    
    % If opts.UseModelInputs is false, the number of variables can change
    if opts.UseModelInputs
        UseControls_i = opts.UseControls{1};
    else
        UseControls_i = opts.UseControls{iCon};
    end
    intOpts.UseControls = UseControls_i;
    inTq = nnz(UseControls_i);
    
    inT = nTk + inTs + inTq;

    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.ObjWeights = opts.ObjWeights(:,iCon);
    tGet = opts.tGet{iCon};
    
    % Integrate
    if opts.continuous(iCon) && opts.complex(iCon)
        sol = integrateObj(m, con(iCon), obj(:,iCon), intOpts);
    elseif opts.complex(iCon)
        sol = integrateSys(m, con(iCon), intOpts);
    elseif opts.continuous(iCon)
        sol = integrateObjSelect(m, con(iCon), obj(:,iCon), tGet, intOpts);
    else
        sol = integrateSysSelect(m, con(iCon), tGet, intOpts);
    end
    
    % Add fields for prior objectives
    sol.UseParams = opts.UseParams;
    sol.UseSeeds = UseSeeds_i;
    sol.UseControls = UseControls_i;
    
    % Extract continuous term
    if opts.continuous(iCon)
        contG = sol.y(nx+1,end);
    else
        contG = 0;
    end
    
    % Compute discrete term
    discG = 0;
    for iObj = 1:nObj
        iDiscG = obj(iObj,iCon).G(sol);
        discG = discG + opts.ObjWeights(iObj,iCon) * iDiscG;
    end
    
    % Add to cumulative goal value
    G = G + contG + discG;
    
    if verboseAll; fprintf('iCon = %d\tG = %g\tTime = %0.2f\n', iCon, contG + discG, toc); end
end

if verbose; fprintf('Summary: G = %g\n', G); end
