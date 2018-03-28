function obs = observationEvents(tF, events, max_events)
%observationEvents Store the simulation at events
%
%   obs = observationEvents(tF, events, max_events)
%
%   This observation scheme monitors a list of events. These events are recorded
%   until either a maximum number of events or a maximum time is encountered.
%
%   Inputs
%   tF: [ nonnegative scalar ]
%       The time at which the simulation will be stopped regardless of how manyu
%       events have fired. (This prevents infinitely long simulations.)
%   events: [ Event vector ne ]
%       The events that will be recorded.
%   max_events: [ nonnegative integer vector ne | nonnegative integer | {inf} ]
%       The maximum number of events to record for each event. If scalar, the
%       max will be broadcast to all events. If infinite (the default), events
%       will be recorded until tF is reached.
%
%   Outputs
%   obs: [ observation struct scalar ]
%       A KroneckerBio observation scheme structure
%
%   For the meanings of the fields of obs see "help observationZero"

% (c) 2018 David R Hagen and Bruce Tidor
% This work is released under the MIT license.

if nargin < 3
    max_events = [];
end

assert(isnumeric(max_events) && (isscalar(max_events) || isempty(max_events) || numel(events) == numel(max_events)), ...
    'KroneckerBio:observationEvents:max_events', ...
    'max_events must be numeric and either a scalar or have same number of elements as events')

if isempty(max_events)
    max_events = inf(size(events));
elseif isscalar(max_events)
    max_events = zeros(size(events)) + max_events;
end

events = vec(events);
max_events = vec(max_events);

ne = numel(events);

obs.Complex = false;

obs.tF = tF;
obs.ne = ne;
obs.Events = @event_stitcher;
obs.IsFinished = @is_finished;

obs.Simulation = @simulation;
obs.Sensitivity = @sensitivity;
obs.Curvature = @curvature;

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        sim.Type = 'Simulation.System.Events';
        sim.Name = int.Name;
        sim.int  = int;
        
        keep_events = mark_events(ne, max_events, int.ie);
        sim.ie   = int.ie(keep_events);
        sim.te   = int.te(keep_events);
        sim.xe   = int.xe(:,keep_events);
        sim.ue   = int.ue(:,keep_events);
        sim.ye   = int.ye(:,keep_events);
    end

    function sim = sensitivity(int)
        sim = simulation(int);
        sim.Type = 'Simulation.Sensitivity.Events';

        keep_events = mark_events(ne, max_events, int.ie);
        sim.dxedT = int.dxedT(:,keep_events);
        sim.duedT = int.duedT(:,keep_events);
        sim.dyedT = int.dyedT(:,keep_events);
    end

    function sim = curvature(int)
        sim = sensitivity(int);
        sim.Type = 'Simulation.Curvature.Events';

        keep_events = mark_events(ne, max_events, int.ie);
        sim.d2xedT2 = int.d2xedT2(:,keep_events);
        sim.d2uedT2 = int.d2uedT2(:,keep_events);
        sim.d2yedT2 = int.d2yedT2(:,keep_events);
    end

    function [value, isterminal, direction] = event_stitcher(t,y)
        value = zeros(ne,1);
        isterminal = max_events ~= inf;
        direction = zeros(ne,1);
        for ie = 1:ne
            value(ie) = events(ie).e(t,y);
            direction(ie) = events(ie).direction;
        end
    end

    function finished = is_finished(sol)
        finished = true;
        for i_event = 1:ne
            if nnz(sol.ie == i_event) < max_events(i_event)
                finished = false;
                break
            end
        end
    end
end

function keep_events = mark_events(ne, max_events, ies)
% Remove events if too many were recorded
keep_events = true(size(ies));
event_counts = zeros(ne,1);
for i = 1:numel(keep_events)
    ie = ies(i);
    event_counts(ie) = event_counts(ie) + 1;
    if event_counts(ie) > max_events(ie)
        keep_events(i) = false;
    end
end
end
