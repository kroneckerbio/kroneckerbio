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
sim.Type = 'Simulation.System.Zero';
sim.Name = int.Name;
sim.int = int;
end

function sim = empty_sensitivity(int)
sim.Type = 'Simulation.Sensitivity.Zero';
sim.Name = int.Name;
sim.int = int;
end

function sim = empty_curvature(int)
sim.Type = 'Simulation.Curvature.Zero';
sim.Name = int.Name;
sim.int = int;
end
