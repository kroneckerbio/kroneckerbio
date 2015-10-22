function int = integrateSensSimpFinite(m, con, tF, eve, fin, t_get, opts)
% Store starting parameter sets
T0 = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, {opts.UseInputControls}, {opts.UseDoseControls});

% Constants
nx = m.nx;
nu = m.nu;
ny = m.ny;
nT = numel(T0);

% Initial simulation
int0 = integrateSysSimp(m, con, tF, eve, fin, t_get, opts);
if numel(int0) > 1; error('Not implemented yet'); end

% Need an observation that watches all discrete times
nt_discrete = numel(int0.t);
nt_events = numel(int0.te);
nt = nt_discrete + nt_events;

t_all = [int0.t, int0.te];
x_all = [int0.x, int0.xe];
u_all = [int0.u, int0.ue];
y_all = [int0.y, int0.ye];

% System results
int.Type = 'Integration.Sensitivity.Simple';
int.Name = [m.Name ' in ' con.Name];

int.x_names = vec({m.States.Name});
int.u_names = vec({m.Inputs.Name});
int.y_names = vec({m.Outputs.Name});
int.k_names = vec({m.Parameters.Name});
int.s_names = vec({m.Seeds.Name});

int.nx = nx;
int.ny = m.ny;
int.nu = m.nu;
int.nk = m.nk;
int.ns = m.ns;
int.nq = con.nq;
int.nh = con.nh;
int.k = m.k;
int.s = con.s;
int.q = con.q;
int.h = con.h;

int.nT = nT;
int.Normalized = opts.Normalized;
int.UseParams = opts.UseParams;
int.UseSeeds = opts.UseSeeds;
int.UseInputControls = opts.UseInputControls;
int.UseDoseControls = opts.UseDoseControls;

int.t = int0.t;
int.x = int0.x;
int.u = int0.u;
int.y = int0.y;

% Change each parameter a bit and resimulate
dxdT_all = zeros(nx*nT,nt);
dudT_all = zeros(nu*nT,nt);
dydT_all = zeros(ny*nT,nt);
for iT = 1:nT
    % Set baseline parameters
    T_i = T0(iT);
    T_up = T0;
    
    % Change current parameter by finite amount
    step_size = 1e-8;
    if opts.Normalized
        norm_factor = T_i;
    else
        norm_factor = 1;
    end
    if opts.ComplexStep
        imag_factor = 1i;
    else
        imag_factor = 1;
    end
    diff = step_size * norm_factor * imag_factor;
    
    % Simulate difference
    T_up(iT) = T_up(iT) + diff;
    [m_up, con_up] = updateAll(m, con, T_up, opts.UseParams, opts.UseSeeds, {opts.UseInputControls}, {opts.UseDoseControls});
    int_up = integrateSysSimp(m_up, con_up, tF, eve, fin, t_all, opts);
    
    % Difference
    ind_x_start = (iT-1)*nx+1;
    ind_x_end = (iT-1)*nx+nx;
    
    ind_u_start = (iT-1)*nu+1;
    ind_u_end = (iT-1)*nu+nu;
    
    ind_y_start = (iT-1)*ny+1;
    ind_y_end = (iT-1)*ny+ny;
    
    if opts.ComplexStep
        dxdT_all(ind_x_start:ind_x_end,:) = imag(int_up.x) ./ step_size;
        dudT_all(ind_u_start:ind_u_end,:) = imag(int_up.u) ./ step_size;
        dydT_all(ind_y_start:ind_y_end,:) = imag(int_up.y) ./ step_size;
    else
        dxdT_all(ind_x_start:ind_x_end,:) = (int_up.x - x_all) ./ step_size;
        dudT_all(ind_u_start:ind_u_end,:) = (int_up.u - u_all) ./ step_size;
        dydT_all(ind_y_start:ind_y_end,:) = (int_up.y - y_all) ./ step_size;
    end
end

int.dxdT = dxdT_all(:,1:nt_discrete);
int.dudT = dudT_all(:,1:nt_discrete);
int.dydT = dydT_all(:,1:nt_discrete);

int.ie = int0.ie;
int.te = int0.te;
int.xe = int0.xe;
int.ue = int0.ue;
int.ye = int0.ye;

int.dxedT = dxdT_all(:,nt_discrete+1:end);
int.duedT = dudT_all(:,nt_discrete+1:end);
int.dyedT = dydT_all(:,nt_discrete+1:end);

int.sol = int0.sol;
