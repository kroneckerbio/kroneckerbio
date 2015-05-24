function obs = observationZero(dims)

if nargin < 1
    dims = 1;
end

dims = row(dims);

if isscalar(dims)
    dims = [dims, 1];
end

obs.Type = 'Observation';
obs.Name = 'UnamedObservation';

obs.Complex = false;

obs.tF = 0;
obs.DiscreteTimes = zeros(1,0);

obs.ne = 0;
obs.Events = @empty_event;
obs.IsFinished = @(sol)true;

obs.Simulation = @empty_simulation;
obs.Sensitivity = @empty_sensitivity;
obs.Curvature = @empty_curvature;

obs.Objective = @(int)objectiveZero();

% Return the observation at the requested size
obs = repmat(obs, dims);
end

function [value, is_terminal, direction] = empty_event(t, y)
value = zeros(0,1);
is_terminal = false(0,1);
direction = zeros(0,1);
end

function sim = empty_simulation(int)
nx = int.nx;
nu = int.nu;
ny = int.ny;

sim.Type = 'Simulation.System.Zero';
sim.Name = int.Name;

sim.t = zeros(1,0);
sim.x = zeros(nx,0);
sim.u = zeros(nu,0);
sim.y = zeros(ny,0);

sim.ie = zeros(1,0);
sim.xe = zeros(nx,0);
sim.ue = zeros(nu,0);
sim.ye = zeros(ny,0);

sim.int = int;
end

function sim = empty_sensitivity(int)
nx = int.nx;
nu = int.nu;
ny = int.ny;
nT = int.nT;

sim = empty_simulation(int);

sim.Type = 'Simulation.Sensitivity.Zero';

sim.dxdT = zeros(nx*nT,0);
sim.dudT = zeros(nu*nT,0);
sim.dydT = zeros(ny*nT,0);

sim.dxedT = zeros(nx*nT,0);
sim.duedT = zeros(nu*nT,0);
sim.dyedT = zeros(ny*nT,0);
end

function sim = empty_curvature(int)
nx = int.nx;
nu = int.nu;
ny = int.ny;
nT = int.nT;

sim = empty_sensitivity(int);

sim.Type = 'Simulation.Curvature.Zero';

sim.d2xdT2 = zeros(nx*nT*nT,0);
sim.d2udT2 = zeros(nu*nT*nT,0);
sim.d2ydT2 = zeros(ny*nT*nT,0);

sim.d2xedT2 = zeros(nx*nT*nT,0);
sim.d2uedT2 = zeros(nu*nT*nT,0);
sim.d2yedT2 = zeros(ny*nT*nT,0);
end
