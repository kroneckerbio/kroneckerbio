function ints = integrateAllSens(m, con, obs, opts, use_finite)
if nargin < 5
    use_finite = false;
end

n_obs = size(obs,1);

nes = vec([obs.ne]);
e_start = cumsum(nes) - nes + 1; % First index of an obs's events
e_end = cumsum(nes); % Last index of an obs's events

[tF, eve, fin, t_get] = collectObservations(m, con, obs);

if any([obs.Complex])
    if ~use_finite
        ints = integrateSensComp(m, con, tF, eve, fin, opts);
    else
        ints = integrateSensCompFinite(m, con, tF, eve, fin, opts);
    end
    
    % Distribute times for each observation
    ints = repmat(ints, n_obs,1);
    for i_obs = 1:n_obs
        if obs(i_obs).Complex
            % Only reveal time points in range of observation
            % Note: deval will still not throw an error outside this range
            ints(i_obs).t = [ints(i_obs).t(ints(i_obs).t < obs(i_obs).tF), obs(i_obs).tF];
        else
            % Evaluate all requested time points
            ints(i_obs).t = obs(i_obs).DiscreteTimes;
            ints(i_obs).x = ints(i_obs).x(ints(i_obs).t);
            ints(i_obs).u = ints(i_obs).u(ints(i_obs).t);
            ints(i_obs).y = ints(i_obs).y(ints(i_obs).t);
            
            ints(i_obs).dxdT = ints(i_obs).dxdT(ints(i_obs).t);
            ints(i_obs).dudT = ints(i_obs).dudT(ints(i_obs).t);
            ints(i_obs).dydT = ints(i_obs).dydT(ints(i_obs).t);
        end
    end
else
    if ~use_finite
        ints = integrateSensSimp(m, con, tF, eve, fin, t_get, opts);
    else
        ints = integrateSensSimpFinite(m, con, tF, eve, fin, t_get, opts);
    end
    
    % Distribute times for each observation
    ints = repmat(ints, n_obs,1);
    for i_obs = 1:n_obs
        % Evaluate all requested time points
        if numel(ints(i_obs).t) ~= numel(obs(i_obs).DiscreteTimes) % Saves time and memory if solution is large
            % lookupmember returns inds_i such that ints(i_obs).t(inds_i) == obs(i_obs).DiscreteTimes
            inds_i = lookupmember(obs(i_obs).DiscreteTimes, ints(i_obs).t);
            inds_i = inds_i(inds_i ~= 0);
            ints(i_obs).t = obs(i_obs).DiscreteTimes;
            ints(i_obs).x = ints(i_obs).x(:,inds_i);
            ints(i_obs).u = ints(i_obs).u(:,inds_i);
            ints(i_obs).y = ints(i_obs).y(:,inds_i);
            
            ints(i_obs).dxdT = ints(i_obs).dxdT(:,inds_i);
            ints(i_obs).dudT = ints(i_obs).dudT(:,inds_i);
            ints(i_obs).dydT = ints(i_obs).dydT(:,inds_i);
        end
    end
end

% Distribute events for each observation
for i_obs = 1:n_obs
    current_events = ints(i_obs).ie >= e_start(i_obs) & ints(i_obs).ie <= e_end(i_obs);
    
    % Change ie to start at 1 for each observation
    ints(i_obs).ie = row(ints(i_obs).ie(current_events));
    ints(i_obs).ie = ints(i_obs).ie - e_start(i_obs) + 1;
    
    ints(i_obs).te = row(ints(i_obs).te(current_events));
    ints(i_obs).xe = ints(i_obs).xe(:,current_events);
    ints(i_obs).ue = ints(i_obs).ue(:,current_events);
    ints(i_obs).ye = ints(i_obs).ye(:,current_events);
    
    ints(i_obs).dxedT = ints(i_obs).dxedT(:,current_events);
    ints(i_obs).duedT = ints(i_obs).duedT(:,current_events);
    ints(i_obs).dyedT = ints(i_obs).dyedT(:,current_events);
end
