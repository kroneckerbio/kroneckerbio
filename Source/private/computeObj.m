function G = computeObj(m, con, obj, opts)
% G = computeObj(m, con, obj, opts)
% This function computes the total objective function value for a vector of
% con and a matrix of obj.
verbose = logical(opts.Verbose);
verbose_all = max(verbose-1,0);

% Constants
nx = m.nx;
n_con = numel(con);
n_obj = size(obj,1);

% Initialize variables
G = 0;

if verbose; disp('Integrating forward...'); end
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
    ints = integrateAllSys(m, con(i_con), obj(:,i_con), opts_i);
    
    % Add fields for prior objectives
    [ints.UseParams] = deal(opts.UseParams);
    [ints.UseSeeds] = deal(UseSeeds_i);
    [ints.UseInputControls] = deal(UseInputControls_i);
    [ints.UseDoseControls] = deal(UseDoseControls_i);
    
    % Extract continuous term
    if opts.continuous(i_con)
        error('Continuous objective functions have not been updated')
        G_cont = ints.y(nx+1,end);
    else
        G_cont = 0;
    end
    
    % Compute discrete term
    G_disc = 0;
    for i_obj = 1:n_obj
        G_disc_i = obj(i_obj,i_con).G(ints(i_obj));
        G_disc = G_disc + opts.ObjWeights(i_obj,i_con) * G_disc_i;
    end
    
    % Add to cumulative goal value
    G = G + G_cont + G_disc;
    
    if verbose_all; fprintf('iCon = %d\tG = %g\tTime = %0.2f\n', i_con, G_cont + G_disc, toc); end
end

if verbose; fprintf('Summary: G = %g\n', G); end
