function int = integrateCurvSimpFinite(m, con, tF, eve, fin, t_get, opts)
% Store starting parameter sets
T0 = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, {opts.UseInputControls}, {opts.UseDoseControls});

% Constants
nx = m.nx;
nu = m.nu;
ny = m.ny;
nT = numel(T0);

% Initial simulation
int0 = integrateSensSimp(m, con, tF, eve, fin, t_get, opts);
if numel(int0) > 1; error('Not implemented yet'); end

% Need an observation that watches all discrete times
nt_discrete = numel(int0.t);
nt_events = numel(int0.te);
nt = nt_discrete + nt_events;

t_all = [int0.t, int0.te];
dxdT_all = [int0.dxdT, int0.dxedT];
dudT_all = [int0.dudT, int0.duedT];
dydT_all = [int0.dydT, int0.dyedT];

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
int.UseParams = opts.UseParams;
int.UseSeeds = opts.UseSeeds;
int.UseInputControls = opts.UseInputControls;
int.UseDoseControls = opts.UseDoseControls;

int.t = int0.t;
int.x = int0.x;
int.u = int0.u;
int.y = int0.y;

int.dxdT = int0.dxdT;
int.dudT = int0.dudT;
int.dydT = int0.dydT;

% Change each parameter a bit and resimulate
d2xdT2_all = zeros(nx*nT*nT,nt);
d2udT2_all = zeros(nu*nT*nT,nt);
d2ydT2_all = zeros(ny*nT*nT,nt);
for iT = 1:nT
    % Set baseline parameters
    Ti = T0(iT);
    T_up = T0;
    
    % Change current parameter by finite amount
    if opts.Normalized
        diff = Ti * 1e-8;
    else
        diff = 1e-8;
    end
    
    % Simulate difference
    T_up(iT) = T_up(iT) + diff;
    [m_up, con_up] = updateAll(m, con, T_up, opts.UseParams, opts.UseSeeds, {opts.UseInputControls}, {opts.UseDoseControls});
    int_up = integrateSensSimp(m_up, con_up, tF, eve, fin, t_all, opts);
    
    % Difference
    ind_xT_start = (iT-1)*nx*nT+1;
    ind_xT_end = (iT-1)*nx*nT+nx*nT;
    d2xdT2_all(ind_xT_start:ind_xT_end,:) = (int_up.dxdT - dxdT_all) ./ diff;

    ind_uT_start = (iT-1)*nu*nT+1;
    ind_uT_end = (iT-1)*nu*nT+nu*nT;
    d2udT2_all(ind_uT_start:ind_uT_end,:) = (int_up.dudT - dudT_all) ./ diff;
    
    ind_yT_start = (iT-1)*ny*nT+1;
    ind_yT_end = (iT-1)*ny*nT+ny*nT;
    d2ydT2_all(ind_yT_start:ind_yT_end,:) = (int_up.dydT - dydT_all) ./ diff;
end

int.d2xdT2 = d2xdT2_all(:,1:nt_discrete);
int.d2udT2 = d2udT2_all(:,1:nt_discrete);
int.d2ydT2 = d2ydT2_all(:,1:nt_discrete);

int.ie = int0.ie;
int.te = int0.te;
int.xe = int0.xe;
int.ue = int0.ue;
int.ye = int0.ye;

int.dxedT = int0.dxedT;
int.duedT = int0.duedT;
int.dyedT = int0.dyedT;

int.d2xedT2 = d2xdT2_all(:,nt_discrete+1:end);
int.d2uedT2 = d2udT2_all(:,nt_discrete+1:end);
int.d2yedT2 = d2ydT2_all(:,nt_discrete+1:end);

int.sol = int0.sol;
