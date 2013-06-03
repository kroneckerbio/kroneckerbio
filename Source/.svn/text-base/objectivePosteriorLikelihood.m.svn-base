function this = objectivePosteriorLikelihood(m, con, obj, opts)
% Clean-up inputs
if nargin < 4
    opts = [];
end

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseICs         = [];
defaultOpts.UseControls    = [];

defaultOpts.ObjWeights     = ones(size(obj));

defaultOpts.Normalized     = true;
defaultOpts.UseAdjoint     = true;

opts = mergestruct(defaultOpts, opts);

opts.Verbose = max(opts.Verbose-1,0);
verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nk = m.nk;
nCon = numel(con);
nObj = size(obj,1);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseICs is a logical matrix
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseControls nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon, m.nq, cat(1,con.nq));

nT = nTk + nTx + nTq;

% Refresh conditions and objectives
con = refreshCon(m, con);
obj = refreshObj(m, con, obj, opts.UseParams, opts.UseICs, opts.UseControls);

% Construct starting variable parameter set
T = collectActiveParameters(m, con, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

discreteTimes = unique(cat(1, 0, opts.tGet{:}));

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 3, opts.continuous, nx, nCon, opts.UseAdjoint, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);

% This objective function is not compatible with the adjoint method
assert(~opts.UseAdjoint, 'KroneckerBio:objectivePosteriorLikelihood:UseAdjoint', 'The Posterior Likelihood objective function is not compatible with the adjoint gradient method.')

% Objective structure
this.Type = 'Objective.Information';
this.Name = 'UnamedObjective';

% Objective function control parameters
this.Continuous = false; % Continuous g here is 0
this.Complex = false; % This objective function does not need it
this.Linked = 0; % This objective function requires all experimental conditions, TODO
this.DiscreteTimes = discreteTimes; % Times at which discrete measurements are taken

% Discrete objective function
this.G       = @G;
this.dGdx    = @dGdx;
this.dGdk    = @dGdk;
this.d2Gdx2  = @d2Gdx2;
this.d2Gdk2  = @d2Gdk2;
this.d2Gdkdx = @d2Gdkdx;
this.d2Gdxdk = @d2Gdxdk;

% Update
this.Update = @Update;

    function [val disc] = G(sol)
        % Compute information matrix for all objectives
        F = zeros(nT,nT);
        intOpts = opts;
        for iCon = 1:nCon
            if size(sol.y, 1) == nx+nx*nT || size(sol.y, 1) == nx+1+nx*nT+nT
                % Correct size
            else%if size(sol.y, 1) == nx || size(sol.y, 1) == nx+1
                % Modify opts structure
                intOpts.AbsTol = opts.AbsTol{iCon};
                intOpts.ObjWeights = opts.ObjWeights(:,iCon);
                tGet = opts.tGet{iCon};
                
                % Sensitivity is required for this objective function
                if opts.continuous(iCon) && opts.complex(iCon)
                    sol = integrateObjSens(m, con(iCon), obj(:,iCon), intOpts);
                elseif opts.complex(iCon)
                    sol = integrateSens(m, con(iCon), intOpts);
                elseif opts.continuous(iCon)
                    sol = integrateObjSensSelect(m, con(iCon), obj(:,iCon), tGet, intOpts);
                else
                    sol = integrateSensSelect(m, con(iCon), tGet, intOpts);
                end
            end
            
            for iObj = 1:nObj
                if opts.Normalized
                    F = F + obj(iObj,iCon).Fn(sol, T);
                else
                    F = F + obj(iObj,iCon).F(sol);
                end
            end
        end
        
        % Compute determinant of information matrix
        lambda = infoeig(F);
        lambda(lambda < 1e-16) = 1e-16;
        val = sum(log(lambda));
        
        disc = discreteTimes;
    end

    function val = dGdx(t, sol)
        % TODO: implement opts.UseICs
        val = zeros(nx,1);
    end

    function val = dGdk(t, sol)
        val = zeros(nk,1);
        if t == 0
            F = zeros(nT,nT);
            dFdT = zeros(nT*nT,nT);
            intOpts = opts;
            
            if verbose; fprintf('Integrating double sensitivities:\n'); end
            for iCon = 1:nCon
                % Modify opts structure
                intOpts.AbsTol = opts.AbsTol{iCon};
                intOpts.ObjWeights = opts.ObjWeights(:,iCon);
                tGet = opts.tGet{iCon};
                
                % Double sensitivity is required for this operation
                if opts.continuous(iCon) && opts.complex(iCon)
                    solDblsens = integrateObjDblsens(m, con(iCon), obj(:,iCon), intOpts);
                elseif opts.complex(iCon)
                    solDblsens = integrateDblsens(m, con(iCon), intOpts);
                elseif opts.continuous(iCon)
                    solDblsens = integrateObjDblsensSelect(m, con(iCon), obj(:,iCon), tGet, intOpts);
                else
                    solDblsens = integrateDblsensSelect(m, con(iCon), tGet, intOpts);
                end
                
                for iObj = 1:nObj
                    if opts.Normalized
                        F = F + obj(iObj,iCon).Fn(sol, T);
                        dFdT = dFdT + obj(iObj,iCon).dFndT(solDblsens, T);
                    else
                        F = F + obj(iObj,iCon).F(sol);
                        dFdT = dFdT + obj(iObj,iCon).dFdT(solDblsens);
                    end
                end
            end
            
            % dGdT = tr(F^-1 * dFdT)
            temp = infoinv(F) * reshape(dFdT, nT,nT*nT);
            
            % Extract diagonal along first and second dimensions
            temp = reshape(temp, nT*nT,nT);
            temp = temp(sub2ind([nT,nT], 1:nT,1:nT),:);
            
            % Compute trace by summing along diagonal elements
            temp = vec(sum(temp, 1));
            
            % Map to UseParams
            val(opts.UseParams) = temp;
        end
    end

    function val = d2Gdx2(t, sol)
        error('You do not want to implement this')
    end

    function val = d2Gdk2(t, sol)
        error('You do not want to implement this')
    end

    function val = d2Gdkdx(t, sol)
        error('You do not want to implement this')
    end

    function val = d2Gdxdk(t, sol)
        error('You do not want to implement this')
    end

    function this = Update(m,con,useParams,useICs,useControls)
        opts.UseParams   = useParams;
        opts.UseICs      = useICs;
        opts.UseControls = useControls;
        this = objectivePosteriorLikelihood(m, con, obj, opts);
    end

end