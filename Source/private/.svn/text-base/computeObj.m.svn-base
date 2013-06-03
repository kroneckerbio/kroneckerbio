function G = computeObj(m, con, obj, opts)
% G = computeObj(m, con, obj, opts)
% This function computes the total objective function value for a vector of
% con and a matrix of obj.
verbose = logical(opts.Verbose);
verboseAll = max(verbose-1,0);

% Constants
nx = m.nx;
nCon = numel(con);
nObj = size(obj,1);

% Initialize variables
G = 0;
intOpts = opts;

if verbose; disp('Integrating forward...'); end
for iCon = 1:nCon
    if verboseAll; tic; end
    
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
