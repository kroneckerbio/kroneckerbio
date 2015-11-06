function obs = observationZero(dims)
%observationZero Default structure for KroneckerBio observation schemes
%
%   obs = observationZero(dims)
%
%   This function returns an observation scheme structure that makes no
%   observations. This file serves as documentation for the observation
%   scheme structure.
%
%   Inputs
%   dims: [ nonnegative integer vector ]
%       The size of the blank objective structure array
%
%   Outputs
%   obs: [ observation scheem structure ]
%       .Type [ 'Observation' ]
%       .Name [ string ]
%           An arbitrary name for the observation scheme
%       .Complex [ logical ]
%           True means that solution is needed at all times. False means
%           that it is only needed at DiscreteTimes.
%       .tF [ nonnegative ]
%           The final time of the simulation needed
%       .DiscreteTimes
%           The particular times at which to save the simulation
%       .ne [ nonegative integer ]
%           The number of event handlers
%       .Events [ function handle ]
%           This is a standard Matlab ODE event handler. It accepts (t,y),
%           where t is a time and y is the output vector at that time. It
%           returns [value, is_terminal, direction], where value is a
%           numeric vector of length ne, is_terminal is a logical vector of
%           length ne, and direction is a vector of length ne of value -1,
%           0, or 1.
%       .IsFinished [ handle @(sol) returns logical ]
%           Each time a terminal event is encountered, the solution up to
%           that point is passed to this function. If this function returns
%           true, then the simulator no longer continues integrating for
%           this observation.
%       .Simulation [ handle @(int) returns sim structure ]
%           Used by SimulateSystem. A function that accepts an basic
%           integration and returns a simulation. Besides a Type field and
%           a Name field, the strucutre is arbitrary.
%       .Sensitivity [ handle @(int) returns sim structure ]
%           Used by SimulateSensitivity. A function that accepts an
%           sensitivity integration and returns a simulation.
%       .Curvature [ handle @(int) returns sim structure ]
%           Used by SimulateCurvature. A function that accepts an
%           curvature integration and returns a simulation.
%       .Lna [ handle @(int) returns sim structure ]
%           Used by SimulateLna. A function that accepts an
%           LNA integration and returns a simulation.
%       .Objective [ handle @(measurements) returns obj structure ]
%           A function that accepts a vector of measurements like
%           sim.measurements and returns an objective function.
%       .F [ handle @(int) returns matrix nT by nT ]
%           Accepts an integration and returns the Fisher information
%           matrix.

% (c) 2015 David R Hagen
% This work is released under the MIT license.

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

obs.Lna = @empty_lna;

obs.Objective = @(int)objectiveZero();

obs.F = @F;

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

function sim = empty_lna(int)
sim.Type = 'Simulation.Lna.Zero';
sim.Name = int.Name;
sim.int = int;
end
