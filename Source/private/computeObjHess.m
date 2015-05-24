function [G, D, H] = computeObjHess(m, con, obj, opts)
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
H = zeros(nT,nT);

if opts.Verbose; fprintf('Integrating curvature:\n'); end
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
    ints = integrateAllCurv(m, con(i_con), obj(:,i_con), opts_i);
    
    % Compute continuous term
    % Extract continuous term
    if opts.continuous(i_con)
        G_cont = ints(1).sol.y(nx+1,end);
        
        dGdTStart = nx+1+nx*inT+1;
        dGdTEnd   = nx+1+nx*inT+inT;
        D_cont = ints(1).sol.y(dGdTStart:dGdTEnd,end);
        
        d2GdT2Start = nx+1+nx*inT+nx*inT*inT+1;
        d2GdT2End   = nx+1+nx*inT+nx*inT*inT+inT;
        H_cont = reshape(ints(1).sol.y(d2GdT2Start:d2GdT2End,end), nT,nT);
    else
        G_cont = 0;
        D_cont = zeros(nT,1);
        H_cont = zeros(nT,nT);
    end
    
    % Compute discrete term
    G_disc = 0;
    D_disc = zeros(inT,1); % T_
    H_disc = zeros(inT,inT); % T_T
    
    discrete_times_all = cell(n_obj,1);
    for i_obj = 1:n_obj
        [G_disc_i, temp] = obj(i_obj,i_con).G(ints(i_obj));
        discrete_times_all{i_obj} = row(unique(temp));
        G_disc = G_disc + opts.ObjWeights(i_obj,i_con) * G_disc_i;
        
        Dk_i = obj(i_obj,i_con).dGdk(ints(i_obj));
        Ds_i = obj(i_obj,i_con).dGds(ints(i_obj));
        Dq_i = obj(i_obj,i_con).dGdq(ints(i_obj));
        Dh_i = obj(i_obj,i_con).dGdh(ints(i_obj));
        D_disc_i = [Dk_i(opts.UseParams); Ds_i(UseSeeds_i); Dq_i(UseInputControls_i); Dh_i(UseDoseControls_i)];

        Hk_i = obj(i_obj,i_con).d2Gdk2(ints(i_obj));
        Hs_i = obj(i_obj,i_con).d2Gds2(ints(i_obj));
        Hq_i = obj(i_obj,i_con).d2Gdq2(ints(i_obj));
        Hh_i = obj(i_obj,i_con).d2Gdh2(ints(i_obj));
        H_disc_i = [Hk_i(opts.UseParams,opts.UseParams), zeros(nTk,nT-nTk);
                    zeros(nTs,nTk), Hs_i(UseSeeds_i,UseSeeds_i), zeros(nTs,nTq+nTh);
                    zeros(nTq,nTk+nTs), Hq_i(UseInputControls_i,UseInputControls_i), zeros(nTq,nTh);
                    zeros(nTh,nT-nTh), Hh_i(UseDoseControls_i,UseDoseControls_i)];
        
        n_disc = numel(discrete_times_all{i_obj});
        for i_disc = 1:n_disc
            ti = discrete_times_all{i_obj}(i_disc);
            if obj(i_obj).Complex
                dydT_i = reshape(ints(i_obj).dydT(ti), nx,nT);
                d2ydT2_i = reshape(ints(i_obj).d2ydT2(ti), nx,nT*nT);
            else
                dydT_i = reshape(ints(i_obj).dydT(:,i_disc), nx,nT);
                d2ydT2_i = reshape(ints(i_obj).d2ydT2(:,i_disc), nx,nT*nT);
            end
            
            dGdy_i = row(obj(i_obj,i_con).dGdy(ti, ints(i_obj)));
            d2Gdy2_i = obj(i_obj,i_con).d2Gdy2(ti, ints(i_obj));
            
            % D = sum(dGdy(t) *{y.y} dydT(t), t) + dGdT
            D_disc_i = D_disc_i + vec(dGdy_i * dydT_i); % _y * y_T -> _T -> T_
            
            % H = (d2Gdy2(t) *{y.y} dydT(t) *{y.y}) dydT(t) + dGdy *{y.y} d2ydT2 + d2GdT2
            H_disc_i = H_disc_i + dydT_i.' * d2Gdy2_i * dydT_i + reshape(dGdy_i * d2ydT2_i, nT,nT); % (y_T -> T_y) * y_y * y_T -> T_T
        end
        
        D_disc = D_disc + opts.ObjWeights(i_obj,i_con) * D_disc_i;
        H_disc = H_disc + opts.ObjWeights(i_obj,i_con) * H_disc_i;
    end
    
    % Sum discrete and continuous terms
    Di = zeros(nT,1);
    Hi = zeros(nT,nT);
    
    % Rate parameters are the same
    inds_k = 1:nTk;

    ind_Ts_start = nTk + sum(sum(opts.UseSeeds(:,1:i_con-1)));
    inds_s = ind_Ts_start+1:ind_Ts_start+inTs;
    inds_si = nTk+1:nTk+inTs;
    
    Tqind = nTk + nTs + sum(cellfun(@sum, opts.UseInputControls(1:i_con-1)));
    inds_q = Tqind+1:Tqind+inTq;
    inds_qi = nTk+inTs+1:nTk+inTs+inTq;
    
    Thind = nTk + nTs + nTq + sum(cellfun(@sum, opts.UseDoseControls(1:i_con-1)));
    inds_h = Thind+1:Thind+inTh;
    inds_hi = nTk+inTs+inTq+1:inT;
    
    inds = [inds_k, inds_s, inds_q, inds_h];
    inds_i = [inds_k, inds_si, inds_qi, inds_hi];
    
    Di(inds) = D_cont(inds_i) + D_disc(inds_i);
    Hi(inds,inds) = H_cont(inds_i,inds_i) + H_disc(inds_i,inds_i);
    
    % Add to cumulative goal value
    G = G + G_cont + G_disc;
    D = D + Di;
    H = H + Hi;

    if verbose_all; fprintf('iCon = %d\t|dGdT| = %g\t||d2GdT2|| = %g\tTime = %0.2f\n', i_con, norm(D_cont + D_disc), det(H_cont + Hdisc), toc); end
end

if opts.Verbose; fprintf('Summary: |dGdT| = %g\t||d2GdT2|| = %g\n', norm(D), det(H)); end
