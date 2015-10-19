function [G, D] = computeObjGrad(m, con, obj, opts)
% Switch to Adjoint method if requested
if opts.UseAdjoint
    [G, D] = computeObjSensAdj(m, con, obj, opts);
    return
end

% Continue with forward method
verbose_all = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
ny = m.ny;
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
    
    % Compute continuous term
    % Extract continuous term
    if opts.continuous(i_con)
        G_cont = ints(1).sol.y(nx+1,end);
        
        dGdTStart = nx+1+nx*inT+1;
        dGdTEnd   = nx+1+nx*inT+inT;
        D_cont = ints(1).sol.y(dGdTStart:dGdTEnd,end);
    else
        G_cont = 0;
        D_cont = zeros(nT,1);
    end
    
    % Compute discrete term
    G_disc = 0;
    D_disc = zeros(inT,1); % T_

    discrete_times_all = cell(n_obj,1);
    for i_obj = 1:n_obj
        [G_disc_i, temp] = obj(i_obj,i_con).G(ints(i_obj));
        discrete_times_all{i_obj} = row(unique(temp));
        G_disc = G_disc + opts.ObjWeights(i_obj,i_con) * G_disc_i;
        
        Dk_i = obj(i_obj,i_con).dGdk(ints(i_obj));
        Ds_i = obj(i_obj,i_con).dGds(ints(i_obj));
        Dq_i = obj(i_obj,i_con).dGdq(ints(i_obj));
        Dh_i = obj(i_obj,i_con).dGdh(ints(i_obj));
        if opts.Normalized
            Dk_i = Dk_i.*ints(i_obj).k;
            Ds_i = Ds_i.*ints(i_obj).s;
            Dq_i = Dq_i.*ints(i_obj).q;
            Dh_i = Dh_i.*ints(i_obj).h;
        end
        D_disc_i = [Dk_i(opts.UseParams); Ds_i(UseSeeds_i); Dq_i(UseInputControls_i); Dh_i(UseDoseControls_i)];

        n_disc = numel(discrete_times_all{i_obj});
        for i_disc = 1:n_disc
            ti = discrete_times_all{i_obj}(i_disc);
            if obj(i_obj).Complex
                dydT_i = reshape(ints(i_obj).dydT(ti), ny, inT);
            else
                dydT_i = reshape(ints(i_obj).dydT(:,i_disc), ny, inT);
            end
            
            dGdy_i = row(obj(i_obj,i_con).dGdy(ti, ints(i_obj)));
            
            % D = sum(dGdy(t) *{y.y} dydT(t), t) + dGdT
            D_disc_i = D_disc_i + vec(dGdy_i * dydT_i); % _y * y_T -> _T -> T_
        end
        
        D_disc = D_disc + opts.ObjWeights(i_obj,i_con) * D_disc_i;
    end
    
    % Sum discrete and continuous terms
    Di = zeros(nT,1);
    
    % Rate parameters are the same
    inds_k = 1:nTk;
    Di(inds_k) = D_cont(inds_k) + D_disc(inds_k);
    
    % s parameters are different
    ind_Ts_start = nTk + sum(sum(opts.UseSeeds(:,1:i_con-1)));
    inds_s = ind_Ts_start+1:ind_Ts_start+inTs;
    inds_si = nTk+1:nTk+inTs;
    Di(inds_s) = D_cont(inds_si) + D_disc(inds_si);
    
    % q parameters are different
    Tqind = nTk + nTs + sum(cellfun(@sum, opts.UseInputControls(1:i_con-1)));
    inds_q = Tqind+1:Tqind+inTq;
    inds_qi = nTk+inTs+1:nTk+inTs+inTq;
    Di(inds_q) = D_cont(inds_qi) + D_disc(inds_qi);
    
    % h parameters are different
    Thind = nTk + nTs + nTq + sum(cellfun(@sum, opts.UseDoseControls(1:i_con-1)));
    inds_h = Thind+1:Thind+inTh;
    inds_hi = nTk+inTs+inTq+1:inT;
    Di(inds_h) = D_cont(inds_hi) + D_disc(inds_hi);
    
    % Add to cumulative goal value
    G = G + G_cont + G_disc;
    D = D + Di;
    
    if verbose_all; fprintf('iCon = %d\t|dGdT| = %g\tTime = %0.2f\n', i_con, norm(D_cont + D_disc), toc); end
end

if opts.Verbose; fprintf('Summary: |dGdT| = %g\n', norm(D)); end
