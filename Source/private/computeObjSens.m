function [G, D] = computeObjSens(m, con, obj, opts)
% Switch to Adjoint method if requested
if opts.UseAdjoint
    [G, D] = computeObjSensAdj(m, con, obj, opts);
    return
end

% Continue with forward method
verbose_all = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nTk = nnz(opts.UseParams);
nTs = nnz(opts.UseSeeds);
nTq = nnz(cat(1,opts.UseInputControls{:}));
nTh = nnz(cat(1,opts.UseDoseControls{:}));
nT  = nTk + nTs + nTq + nTh;
n_con = numel(con);
n_obj = size(obj,1);

% Initialize variables
G = 0;
D = zeros(nT,1);

if opts.Verbose; fprintf('Integrating sensitivities:\n'); end
for i_con = 1:n_con
    if verbose_all; tic; end
    
    % Modify opts structure
    opts_i = opts;
    opts_i.AbsTol = opts.AbsTol{i_con};
    opts_i.ObjWeights = opts.ObjWeights(:,i_con);
    
    UseSeeds_i = opts.UseSeeds(:,i_con);
    opts_i.UseSeeds = UseSeeds_i;
    inTs = nnz(UseSeeds_i);
    
    UseInputControls_i = opts.UseInputControls{i_con};
    opts_i.UseInputControls = UseInputControls_i;
    inTq = nnz(UseInputControls_i);
    
    UseDoseControls_i = opts.UseDoseControls{i_con};
    opts_i.UseDoseControls = UseDoseControls_i;
    inTh = nnz(UseDoseControls_i);
    
    inT = nTk + inTs + inTq + inTh;

    % Integrate
    ints = integrateAllSens(m, con(i_con), obj(:,i_con), opts_i);
    
    % *Compute G*
    % Extract continuous term
    if opts.continuous(i_con)
        G_cont = ints(1).sol.y(nx+1,end);
    else
        G_cont = 0;
    end
    
    % Compute discrete term
    G_disc = 0;
    discrete_times_all = cell(n_obj,1);
    for i_obj = 1:n_obj
        [iDiscG, temp] = obj(i_obj,i_con).G(ints(i_obj));
        discrete_times_all{i_obj} = row(unique(temp));
        G_disc = G_disc + opts.ObjWeights(i_obj,i_con) * iDiscG;
    end
    
    % Add to cumulative goal value
    G = G + G_cont + G_disc;
    
    % *Compute D*
    % Extract continuous term
    if opts.continuous(i_con)
        dGdTStart = nx+1+nx*inT+1;
        dGdTEnd   = nx+1+nx*inT+inT;
        D_cont = ints(1).sol.y(dGdTStart:dGdTEnd,end);
    else
        D_cont = zeros(nT,1);
    end
    
    % Compute discrete term
    D_disc = zeros(inT,1); % T_t
    for i_obj = 1:n_obj
        n_disc = numel(discrete_times_all{i_obj});
        D_disc_i = zeros(nT,1);
        for i_disc = 1:n_disc
            ti = discrete_times_all{i_obj}(i_disc);
            if obj(i_obj).Complex
                dGdx_i = row(obj(i_obj,i_con).dGdx(ti, ints(i_obj)));
                dxdT_i = reshape(ints(i_obj).dxdT(ti), nx,nT);
            else
                dGdx_i = row(obj(i_obj,i_con).dGdx(ti, ints(i_obj)));
                dxdT_i = reshape(ints(i_obj).dxdT(:,i_disc), nx,nT);
            end
            D_disc_i = D_disc_i + vec(dGdx_i * dxdT_i);
        end
        
        D_disc = D_disc + opts.ObjWeights(i_obj,i_con) * D_disc_i;
    end
    
    % Sum discrete and continuous terms
    Di = zeros(nT,1);
    
    % Rate parameters are the same
    Di(1:nTk) = Di(1:nTk) + D_cont(1:nTk) + D_disc(1:nTk);
    
    % s parameters are different
    Tsind = nTk + sum(sum(opts.UseSeeds(:,1:i_con-1)));
    Di(Tsind+1:Tsind+inTs) = D_cont(nTk+1:nTk+inTs) + D_disc(nTk+1:nTk+inTs);
    
    % q parameters are different
    Tqind = nTk + nTs + sum(cellfun(@sum, opts.UseInputControls(1:i_con-1)));
    Di(Tqind+1:Tqind+inTq) = D_cont(nTk+inTs+1:nTk+inTs+inTq) + D_disc(nTk+inTs+1:nTk+inTs+inTq);
    
    % h parameters are different
    Thind = nTk + nTs + nTq + sum(cellfun(@sum, opts.UseDoseControls(1:i_con-1)));
    Di(Thind+1:Thind+inTh) = D_cont(nTk+inTs+inTq+1:end) + D_disc(nTk+inTs+inTq+1:end);
    
    % Update sum
    D = D + Di;
    
    if verbose_all; fprintf('iCon = %d\t|dGdT| = %g\tTime = %0.2f\n', i_con, norm(D_cont + D_disc), toc); end
end

if opts.Verbose; fprintf('Summary: |dGdT| = %g\n', norm(D)); end
