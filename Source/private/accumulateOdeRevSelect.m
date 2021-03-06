function cum_sol = accumulateOdeRevSelect(der, jac, t0, tF, ic, discontinuities, t_get, nonnegative, RelTol, AbsTol, delta, events, t_first_event, n_events)
% function cumSol = accumulateOdeRevSelect(der, jac, t0, tF, ic, u, discontinuities,
% t_get, nonnegative, RelTol, AbsTol, delta, events, t_first_event, n_events)

% All discontinuities must be equal to their limit from right
% All events must be terminal

%% Work-up
% Clean up inputs
if nargin < 14
    n_events = [];
    if nargin < 13
        t_first_event = [];
        if nargin < 12
            events = [];
            if nargin < 11
                delta = [];
            end
        end
    end
end

% Default inputs
if isempty(t_first_event)
    t_first_event = 0;
end
if isempty(n_events)
    n_events = inf;
end

% Set flag to return a differentiable solution or just particular points
t_get = vec(t_get).';

% discontinuities will be a column vector sorted ascending
discontinuities = discontinuities((discontinuities > t0) & (discontinuities < tF));
discontinuities = unique([vec(discontinuities); t0; tF]);

% Performance constants
forward_boost           = 10;   % Number of epsilons to jump at a discontinuity
max_false_events        = 1e3;  % Number at which integrator should give up
min_time_between_events = 1e-3; % Events are invalid if they occur within this timespan
min_species_change      = 1e-6; % Events are invalid if all species change less than this
N = numel(discontinuities);

% ODE options
sim_opts = odeset('Jacobian', jac, 'AbsTol', AbsTol, 'RelTol', RelTol, 'NonNegative', nonnegative);
if ~isempty(events)
    sim_opts.Events = events;
end

%% Initialize variables
% Initialize event counters
ne = 0;             % number of events past t_first_event
ne_all = 0;         % total number of notable events
n_false_events = 0; % total number of false events

% Initialize empty solution
cum_sol.solver = 'ode15s';
cum_sol.extdata = [];
cum_sol.x = t_get;
cum_sol.y = zeros(numel(ic), numel(t_get));
cum_sol.stats = [];
cum_sol.idata = [];

% Apply the appropriate delta for the first discontinuities
if ~isempty(delta)
    ic = ic - delta(discontinuities(N), ic);
end

%% Integrate through each time interval
for k = 1:(N-1)
    % Update the indexes of the time points
    i0 = N-(k-1);
    i1 = N-k;
    
    % Integration time interval between discontinuities
    if i0 == N
        % At final interval evaluate both points
        t_int = [discontinuities(i0), discontinuities(i1)];
    else
        % At intermediate intervals, do not evaluate the forwardmost point
        % Instead, stop a little before going forward or start a little
        % going backward.
        t_int = [discontinuities(i0) - forward_boost*eps(discontinuities(i0)), discontinuities(i1)];
    end
    
    % Initialize variables for tracking events
    t_event_int = t_int;
    t_event_start = t_int(1); % Marks point where to start after an event
    t_event_end = t_int(2); % Does not get updated when events fire
    
    % Integrate until the total number of events is equal to n_events or
    % this time interval is completed
    % Note: t_event_start <= t_event_end if boost goes past the end time
    while t_event_start > t_event_end && ne < n_events && n_false_events < max_false_events
        % Find tGet points for this interval
        it_get = t_get(t_get >= min(t_event_int) & t_get <= max(t_event_int));
        
        % Construct time vector for ode solver
        t_ode = unique([t_event_int, it_get]);
        if numel(t_ode) < 3
            % Include a third element to prevent the solver from saving
            % all timepoints
            t_ode = [t_ode(1), (t_ode(1)+t_ode(2))/2, t_ode(2)];
        end
        t_ode = t_ode(end:-1:1);% Unique will have sorted tOde the wrong way. Reverse order
        
        % Start or restart integration after an event
        if isempty(events)
            [sim_x, sim_y] = ode15sf(der, t_ode, ic, sim_opts);
            sim_sol.x = sim_x';
            sim_sol.y = sim_y';
        else
            [sim_x, sim_y, sim_xe, sim_ye, sim_ie] = ode15sf(der, t_ode, ic, sim_opts);
            sim_sol.x = sim_x';
            sim_sol.y = sim_y';
            sim_sol.xe = sim_xe';
            sim_sol.ye = sim_ye';
            sim_sol.ie = sim_ie';
        end
        
        % Check if integration failed when it shouldn't have
        if ~isempty(events) && isempty(sim_sol.ie) && sim_sol.x(end) ~= t_int(2)
            try evalin('caller', 'm_k = m.k; m_s = m.s; m_q = m.q; save(''odefail.mat'',''m_k'',''m_s'',''m_q'',''con'')'); end
            error('KroneckerBio:accumulateOde:IntegrationFailure', 'Did not integrate through entire interval!');
        end
        
        % Discard all but the last event to prevent repeats
        if isfield(sim_sol, 'ie') && numel(sim_sol.ie) > 1
            sim_sol.ie = sim_sol.ie(end);
            sim_sol.xe = sim_sol.xe(end);
            sim_sol.ye = sim_sol.ye(:,end);
        end
        
        % Handle delta
        if ~isempty(delta)
            sim_sol.y(:,end) = sim_sol.y(:,end) - delta(sim_sol.x(end), sim_sol.y(:,end));
        end
        
        % Update ICs
        ic  = sim_sol.y(:,end);
        
        % Combine new solution with cumulative
        cum_sol = pasteSols(cum_sol, sim_sol);
        
        % Update time interval that remains to be integrated
        t_event_start = sim_sol.x(end) - forward_boost*eps(sim_sol.x(end));
        t_event_int   = [t_event_start, t_event_end];
        
        % If we have not reached tF, then a terminating event was
        % detected. Check to ensure that this is not due to numerical
        % problems, but rather that it constitutes a 'true' event.
        %
        % There are at least two criteria that can invalidate an event:
        % 1 -   It is next to the previous event (i.e. a true event that
        %       causes termination twice in succession at very nearby
        %       timepoints)
        % 2 -   It is due to very small change in a state state variable,
        %       a 'wiggle' with low amplitude.
        
        if sim_sol.x(end) ~= t_int(2)
            % **** IMPORTANT SECTION ****
            % This line makes sure events are not double-counted. This is
            % done by checking if two of the SAME event are within
            % min_time_between_events time units of one another.
            
            % Concatenate the initial conditions to the front
            augmented_event_times = [0, cum_sol.xe];
            augmented_event_amp   = [cum_sol.y(:,1), cum_sol.ye];
            
            % Check for a duplicated event
            is_duplicate = numel(cum_sol.ie) > 1 && ...                 % there has been more than one event
                ~diff(cum_sol.ie(end-1:end)) && ...                     % AND it's the same type as the previous one
                diff(cum_sol.xe(end-1:end) < min_time_between_events);  % AND it's really near the previous one in time
            
            % Check for a wiggle event
            is_wiggle = diff(augmented_event_times(end-1:end)) < min_time_between_events || ...   % this event happened really near the previous one in time
                all(abs(diff(augmented_event_amp(:,end-1:end), [], 2)) < min_species_change);     % AND no species changed appreciably
            
            if is_duplicate || is_wiggle
                % Discard the duplicate event
                cum_sol.xe = cum_sol.xe(1:ne_all);
                cum_sol.ie = cum_sol.ie(1:ne_all);
                cum_sol.ye = cum_sol.ye(:,1:ne_all);
                
                % Record its duplicity
                n_false_events = n_false_events + 1;
            else
                % Note the event
                ne_all = ne_all + 1;
                if sim_sol.x(end) > t_first_event
                    % Note the real event
                    ne = ne + 1;
                end
            end
        end
    end
end

%% Work-down
if N == 1
    % If no integration was performed:
    % t0 == tF == discontinuities
    sim_sol.x = t0; % To prevent NaN checker from crashing
    cum_sol.y(:,t0 == t_get) = repmat(ic, 1,nnz(t0 == t_get)); % Insert initial conditions appropriately
end

% Make non-events solutions compatible with those that expect them
if isempty(events)
    cum_sol.xe = [];
    cum_sol.ye = [];
    cum_sol.ie = [];
end

% Put NaN in columns that were never reached during integration
times_never_reached = cum_sol.x < sim_sol.x(end) | cum_sol.x > tF;
cum_sol.y(:,times_never_reached) = NaN;

% Order solution structure fields (in case cumSol.xe, .ye and .ie are out of order)
cum_sol = orderfields(cum_sol, {'solver', 'extdata', 'x', 'y',  'xe', 'ye', 'ie', 'stats', 'idata'});
end