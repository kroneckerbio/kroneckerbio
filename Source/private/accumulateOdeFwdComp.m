function cum_sol = accumulateOdeFwdComp(der, jac, t0, tF, ic, discontinuities, nonnegative, RelTol, AbsTol, delta, events, is_finished)
% function cumSol = accumulateOdeFwdCom[(der, jac, t0, tF, ic, discontinuities, 
% nonnegative, RelTol, AbsTol, delta, events, is_finished)

% All discontinuities must be equal to their limit from right

%% Work-up
% Clean up inputs
if nargin < 12
    is_finished = [];
    if nargin < 11
        events = [];
        if nargin < 10
            delta = [];
        end
    end
end

% discontinuities will be a column vector sorted ascending
discontinuities = discontinuities((discontinuities > t0) & (discontinuities < tF));
discontinuities = unique([vec(discontinuities); t0; tF]);
N = numel(discontinuities);

% Performance constants
forward_boost = 4;   % Number of epsilons to jump at a discontinuity

% ODE options
sim_opts = odeset('Jacobian', jac, 'AbsTol', AbsTol, 'RelTol', RelTol, 'NonNegative', nonnegative);
if ~isempty(events)
    sim_opts.Events = events;
end

%% Initialize variables
% Initialize empty solution
cum_sol = [];

% Apply the appropriate delta for the first discontinuities
if ~isempty(delta)
    ic = ic + delta(discontinuities(1), ic);
end

%% Integrate through each time interval
finished = false;
for k = 1:(N-1)
    % Update the indexes of the time points
    i0 = k;
    i1 = k+1;
    
    % Integration time interval between discontinuities
    if isinf(discontinuities(i1))
        t_int = [discontinuities(i0), discontinuities(i1)];
    else
        t_int = [discontinuities(i0), discontinuities(i1) - forward_boost*eps(discontinuities(i1))];
    end
    
    % Initialize variables for tracking events
    t_event_int = t_int;
    
    % Integrate until interval is complete or terminal event fires
    while true %dowhile
        % Start or restart integration after an event
        sim_sol = ode15sf(der, t_event_int, ic, sim_opts);
        
        % Check if integration failed when it shouldn't have
        if isempty(events) && sim_sol.x(end) ~= t_int(2)
            try evalin('caller', 'm_k = m.k; m_s = m.s; m_q = m.q; save(''odefail.mat'',''m_k'',''m_s'',''m_q'',''con'')'); end
            error('KroneckerBio:accumulateSol:IntegrationFailure', 'Did not integrate through entire interval!');
        end
        
        % Handle delta
        if ~isempty(delta) && sim_sol.x(end) == t_int(2)
            % Deltas are evaluated with boost removed only if no events triggered
            ic = sim_sol.y(:,end) + delta(discontinuities(i1), sim_sol.y(:,end));
        else
            ic = sim_sol.y(:,end);
        end
        
        % Combine new solution with cumulative
        cum_sol = mergeSols(cum_sol, sim_sol);
        
        
        % Test for ending simulation if stop was caused by event
        if sim_sol.x(end) == t_int(2)
            % Segment complete; continue with next segment.
            break
        elseif ~isempty(is_finished) && is_finished(cum_sol)
            % Terminal event triggered
            finished = true;
            break;
        else
            % Update time interval that remains to be integrated
            t_event_int = [sim_sol.x(end) + forward_boost*eps(sim_sol.x(end)), t_int(2)];
            %continue
        end
    end
    
    % Test if last segment finished simulation
    if finished
        break;
    end
end

%% Work-down
if N == 1
    % If no integration was performed:
    % t0 == tF == discontinuities
    cum_sol.solver = 'ode15s';
    cum_sol.extdata = [];
    cum_sol.x = t0;
    cum_sol.y = ic;
    cum_sol.stats = [];
    cum_sol.idata.kvec = 0;
    cum_sol.idata.dif3d = zeros(N,1,1);
    cum_sol.idata.idxNonNegative = sim_opts.NonNegative;
else
    % Handle final point (mainly because of delta at last point)
    if all(ic == sim_sol.y(:,end))
        % Just assume that we made it to the end
        if ~isinf(tF)
            % Don't replace final point with Inf because this causes
            % interpolation errors
            cum_sol.x(end) = tF;
        end
    else
        % A new point has to be added to the end to handle discontinuities
        % or t=Inf
        sim_sol = struct;
        sim_sol.solver = 'ode15s';
        sim_sol.extdata = [];
        sim_sol.x = tF;
        sim_sol.y = ic;
        sim_sol.stats = [];
        sim_sol.idata.kvec = 0;
        sim_sol.idata.dif3d = zeros(size(cum_sol.idata.dif3d,1),1,1);
        sim_sol.idata.idxNonNegative = sim_opts.NonNegative;
        
        cum_sol = mergeSols(cum_sol, sim_sol);
    end
end

% Make non-events solutions compatible with those that expect them
if isempty(events) || ~isfield(cum_sol, 'ie') || isempty(cum_sol.ie)
    cum_sol.xe = zeros(1,0);
    cum_sol.ye = zeros(size(cum_sol.y,1),0);
    cum_sol.ie = zeros(1,0);
end

% Order solution structure fields (in case cumSol.xe, .ye and .ie are out of order)
cum_sol = orderfields(cum_sol, {'solver', 'extdata', 'x', 'y', 'xe', 'ye', 'ie', 'stats', 'idata'});
end
