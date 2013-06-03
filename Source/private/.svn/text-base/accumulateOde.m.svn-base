function cumSol = accumulateOde(der, jac, t0, tF, ic, u, discontinuities, nonnegative, RelTol, AbsTol, delta, direction, events, tFirstEvent, nEvents, tGet)
% function cumSol = accumulateOde(der, jac, t0, tF, ic, u, discontinuities,
% nonnegative, RelTol, AbsTol, delta, direction, events, tFirstEvent, nEvents, tGet)

% All discontinuities must be equal to their right limit
% All events must be terminal
% Events must only exist in forward integration

%% Work-up
% Clean up inputs
if nargin < 16
    tGet = [];
    if nargin < 15
        nEvents = [];
        if nargin < 14
            tFirstEvent = [];
            if nargin < 13
                events = [];
                if nargin < 12
                    direction = [];
                    if nargin < 11
                        delta = [];
                    end
                end
            end
        end
    end
end

% Default inputs
if isempty(tFirstEvent)
    tFirstEvent = 0;
end
if isempty(nEvents)
    nEvents = inf;
end
if isempty(direction)
    direction = 1;
end

% Set flag to return a differentiable solution or just particular points
tGet = vec(tGet).';
if isempty(tGet)
    fullSol = true;
else
    fullSol = false;
end

% discontinuities will be a column vector sorted ascending
discontinuities = discontinuities((discontinuities > t0) & (discontinuities < tF));
discontinuities = unique([vec(discontinuities); t0; tF]);

% Error checking inputs
assert(abs(direction) == 1, 'KroneckerBio:accumulateSol', 'direction must be 1 (forward) or -1 (reverse)');

% Performance constants
forwardBoost            = 10;   % Number of epsilons to jump at a discontinuity
maxFalseEvents          = 1e3;  % Number at which integrator should give up
minTimeBetweenEvents    = 1e-3; % Events are invalid if they occur within this timespan
minSpeciesChange        = 1e-6; % Events are invalid if all species change less than this
N = length(discontinuities);

% ODE options
simOpts = odeset('Jacobian', jac, 'AbsTol', AbsTol, 'RelTol', RelTol, 'NonNegative', nonnegative);
if ~isempty(events)
    simOpts.Events = events;
end

%% Initialize variables
% Initialize event counters
ne = 0;             % number of events past tFirstEvent
neAll = 0;          % total number of notable events
nFalseEvents = 0;   % total number of false events

% Initialize empty solution
if fullSol
    cumSol = [];
else
    cumSol.solver = 'ode15s';
    cumSol.extdata = [];
    cumSol.x = tGet;
    cumSol.y = zeros(length(ic),length(tGet));
    cumSol.stats = [];
    cumSol.idata = [];
end

% Apply the appropriate delta for the first discontinuities
if ~isempty(delta)
    if direction == 1
        ic = ic + direction*delta(discontinuities(1));
    else%if direction == -1
        ic = ic + direction*delta(discontinuities(N));
    end
end

%% Integrate through each time interval
for k = 1:(N-1)
    % Update the indexes of the time points
    if direction == 1   	% Forward integration
        i0 = k;
        i1 = k+1;
    else%if direction == -1 % Reverse integration
        i0 = N-(k-1);
        i1 = N-k;
    end
    
    % Integration time interval
    if i0 == N || i1 == N
        % At final interval evaluate both points
        tInt = [discontinuities(i0), discontinuities(i1)];
    else
        % At intermediate intervals, do not evaluate the forwardmost point
        % Instead, stop a little before going forward or start a little
        % going backward.
        if direction == 1       % Forward
            tInt = [discontinuities(i0), discontinuities(i1) - forwardBoost*eps(discontinuities(i1))];
        else%if direction == -1 % Reverse
            tInt = [discontinuities(i0) - forwardBoost*eps(discontinuities(i0)), discontinuities(i1)];
        end
    end
    
    % Initialize variables for tracking events
    tEventInt = tInt;
    t0 = tInt(1);
    
    % Integrate until the total number of events is equal to nEvents or
    % this time interval is completed
    while ne < nEvents && t0 ~= tEventInt(2) && nFalseEvents < maxFalseEvents
        if fullSol
            % We want the full solution in the interval tEventInt
            tOde = tEventInt;
        else
            % Find tGet points for this interval
            if tEventInt(1) == discontinuities(end) || tEventInt(2) == discontinuities(end)
                % Find points matching the end
                itGet = tGet(tGet >= min(tEventInt) & tGet <= max(tEventInt));
            else
                % The end point is covered by another interval
                itGet = tGet(tGet >= min(tEventInt) & tGet < max(tEventInt));
            end
            
            % Construct time vector for ode solver
            tOde = unique([tEventInt, itGet]);
            if numel(tOde) < 3
                % Include a third element to prevent the solver from saving
                % all timepoints
                tOde = [tOde(1), (tOde(1)+tOde(2))/2, tOde(2)];
            end
            if direction == -1
                % Unique will have sorted tOde the wrong way. Reverse order
                tOde = tOde(end:-1:1);
            end
        end
        
        % Start or restart integration after an event
        if fullSol
            % Ask for single output
            simSol = ode15sf(der, tOde, ic, simOpts, u);
        else
            % Multiple outputs are required
            if isempty(events)
                [simx simy] = ode15sf(der, tOde, ic, simOpts, u);
                simSol.x = simx';
                simSol.y = simy';
            else
                [simx simy simxe simye simie] = ode15sf(der, tOde, ic, simOpts, u);
                simSol.x = simx';
                simSol.y = simy';
                simSol.xe = simxe';
                simSol.ye = simye';
                simSol.ie = simie';
            end
        end
        
        % Check if integration failed when it shouldn't have
        if isempty(events) && simSol.x(end) ~= tInt(2)
            try evalin('caller', 'save(''odefail.mat'',''m'',''con'')'); end
            error('KroneckerBio:accumulateSol:IntegrationFailure', 'Did not integrate through entire interval!');
        end
        
        % Discard all but the last event to prevent repeats
        if isfield(simSol, 'ie') && length(simSol.ie) > 1
            simSol.ie = simSol.ie(end);
            simSol.xe = simSol.xe(end);
            simSol.ye = simSol.ye(:,end);
        end
        
        % Combine new solution with cumulative
        if fullSol
            cumSol = mergeSols(cumSol, simSol);
        else
            cumSol = pasteSols(cumSol, simSol);
        end
        
        % Handle any discrete times or event dependent times
        if ~isempty(delta)
            ic = simSol.y(:, end) + direction*delta(simSol.x(end));
        else
            ic  = simSol.y(:, end);
        end
        
        % Update time interval that remains to be integrated
        t0         = simSol.x(end);
        tEventInt  = [t0 + forwardBoost*direction*eps(t0), tInt(2)];
        
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
        
        if simSol.x(end) ~= tInt(2)
            % **** IMPORTANT SECTION ****
            % This line makes sure events are not double-counted. This is done
            % by checking if two of the SAME event are within 1e-5 time units
            % of one another.
            
            % Concatenate the initial conditions to the front
            augmentedEventTimes = [0 cumSol.xe];
            augmentedEventAmp   = [cumSol.y(:,1) cumSol.ye];
            
            % Check for a duplicated event
            isDuplicate = length(cumSol.ie) > 1 && ...              % there has been more than one event
                ~diff(cumSol.ie(end-1:end)) && ...                  % AND it's the same type as the previous one
                diff(cumSol.xe(end-1:end) < minTimeBetweenEvents);  % AND it's really near the previous one in time
            
            % Check for a wiggle event
            isWiggle = diff(augmentedEventTimes(end-1:end)) < minTimeBetweenEvents || ...     % this event happened really near the previous one in time
                all(abs(diff(augmentedEventAmp(:,end-1:end), [], 2)) < minSpeciesChange);     % AND no species changed appreciably
            
            if isDuplicate || isWiggle
                % Discard the duplicate event
                cumSol.xe = cumSol.xe(1:neAll);
                cumSol.ie = cumSol.ie(1:neAll);
                cumSol.ye = cumSol.ye(:,1:neAll);
                
                % Record its duplicity
                nFalseEvents = nFalseEvents + 1;
            else
                % Note the event
                neAll = neAll + 1;
                if simSol.x(end) > tFirstEvent
                    % Note the real event
                    ne = ne + 1;
                end
            end
        end
    end
end

%% Work-down
% If no integration was performed return an appropriate structure
if N == 1
    if fullSol
        cumSol.solver = 'ode15s';
        cumSol.extdata = [];
        cumSol.x = discontinuities;
        cumSol.y = ic;
        cumSol.stats = [];
        cumSol.idata = [];
    else
        simSol.x = discontinuities;
        cumSol.y(:,discontinuities == tGet) = repmat(ic, 1,nnz(discontinuities == tGet));
    end
end

% Add on any changes to initial conditions at the final timepoint
if ~isempty(delta)
    if fullSol
        cumSol.y(:, end) = cumSol.y(:, end) + direction*delta(cumSol.x(end));
    else
        % Find the final indexes and add to the solution at those points
        finalInds = find(simSol.x(end) == cumSol.x);
        cumSol.y(:,finalInds) = cumSol.y(:,finalInds) + repmat(direction*delta(simSol.x(end)), 1, length(finalInds));
    end
end

% Make non-events solutions compatible with those that expect them
if isempty(events)
    cumSol.xe = [];
    cumSol.ye = [];
    cumSol.ie = [];
end

% Put NaN in columns that were never reached during integration
if ~fullSol
    if direction == 1
        cumSol.y(:,cumSol.x > simSol.x(end)) = NaN;
    elseif direction == -1
        cumSol.y(:,cumSol.x < simSol.x(end)) = NaN;
    end
end

% Order solution structure fields (in case cumSol.xe, .ye and .ie are out of order)
cumSol = orderfields(cumSol, {'solver', 'extdata', 'x', 'y',  'xe', 'ye', 'ie', 'stats', 'idata'});
end