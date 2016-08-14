function [G,D] = computeError(m, con, obj, opts)
% G = computeObj(m, con, obj, opts)
% This function computes the total objective function value for a vector of
% con and a matrix of obj.
verbose = logical(opts.Verbose);
verbose_all = max(verbose-1,0);
calcJac = nargout > 1;

% Constants
n_con = numel(con);
n_obj = size(obj,1);

G = zeros(100,1);
if calcJac
    nTk = nnz(opts.UseParams);
    nTs = nnz(opts.UseSeeds);
    nTq = nnz(cat(1,opts.UseInputControls{:}));
    nTh = nnz(cat(1,opts.UseDoseControls{:}));
    nT  = nTk + nTs + nTq + nTh;
    D = zeros(100,nT);
end

if verbose; disp('Integrating forward...'); end
count = 1;
for i_con = 1:n_con
    if verbose_all; tic; end
    
    % Modify opts structure
    opts_i = opts;
    
    UseSeeds_i = opts.UseSeeds(:,i_con);
    opts_i.UseSeeds = UseSeeds_i;
    
    UseInputControls_i = opts.UseInputControls{i_con};
    opts_i.UseInputControls = UseInputControls_i;
    
    UseDoseControls_i = opts.UseDoseControls{i_con};
    opts_i.UseDoseControls = UseDoseControls_i;
    
    opts_i.AbsTol = opts.AbsTol{i_con};
    opts_i.ObjWeights = opts.ObjWeights(:,i_con);
    
    % Integrate
    if nargout == 1
        ints = integrateAllSys(m, con(i_con), obj(:,i_con), opts_i);
    elseif nargout == 2
        ints = integrateAllSens(m, con(i_con), obj(:,i_con), opts_i);
    end
    
    % Add fields for prior objectives
    [ints.UseParams] = deal(opts.UseParams);
    [ints.UseSeeds] = deal(UseSeeds_i);
    [ints.UseInputControls] = deal(UseInputControls_i);
    [ints.UseDoseControls] = deal(UseDoseControls_i);
    
    if calcJac
        [ints.nT] = deal(nT);
    end
    
    count_start_icon = count;
    for i_obj = 1:n_obj
        G_i = obj(i_obj,i_con).err(ints(i_obj));
        findex = count+numel(G_i)-1;
        addRows = findex > numel(G);
        if addRows
            factor2 = ceil(log2(findex/numel(G)));
            addsize = numel(G)*(2^factor2-1);
            G = [G; zeros(addsize,1)];
        end
        G(count:findex) = G(count:findex) + opts.ObjWeights(i_obj,i_con) .* G_i;
        
        if calcJac
            D_i = obj(i_obj,i_con).derrdT(ints(i_obj));
            if addRows
                D = [D; zeros(addsize,nT)];
            end
            D(count:findex,:) = D(count:findex,:) + opts.ObjWeights(i_obj,i_con) .* D_i;
        end
        
        count = count + numel(G_i);
    end
    
    if verbose_all; fprintf('iCon = %d\tResidual = %g\tTime = %0.2f\n', i_con, sum(G(count_start_icon:findex).^2), toc); end
end

G = G(1:findex);
if calcJac
    D = D(1:findex,:);
end

if verbose; fprintf('Summary: Residual = %g\n', sum(G.^2)); end
