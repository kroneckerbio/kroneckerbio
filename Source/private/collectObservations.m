function [tF, eve, fin, t_get] = collectObservations(m, con, obs)
nx = m.nx;
u = con.u;
y = m.y;

n_obs = numel(obs);
nes = vec([obs.ne]);
n_eve = sum(nes);

tF = max([obs.tF]);
eve = @events_combined;
fin = @is_finished_combined;
t_get = unique([obs.DiscreteTimes]);

    function [val, is_terminal, direction] = events_combined(t, joint)
        x = joint(1:nx);
        
        val = zeros(n_eve,1);
        is_terminal = zeros(n_eve,1);
        direction = zeros(n_eve,1);
        
        u_t = u(t);
        y_t = y(t,x,u_t);
        
        i_eve = 0;
        for iobs = 1:n_obs
            i_start = i_eve + 1;
            i_end = i_eve + nes(iobs);
            i_int = i_start:i_end;
            [val(i_int), is_terminal(i_int), direction(i_int)] = obs(iobs).Events(t,y_t);
        end
    end

    % Is finished
    function finished = is_finished_combined(sol)
        finished = true;
        used_events = 0;
        for iobs = 1:n_obs
            % In the full list of events, which ie values pertain to this obs
            obs_event_indexes = (1:obs(iobs).ne) + used_events;
            
            % Slice out just the events pertaining to this obs
            keep_indexes = ismember(sol.ie, obs_event_indexes);
            sol_i = sol;
            sol_i.ie = sol_i.ie(keep_indexes);
            sol_i.xe = sol_i.xe(keep_indexes);
            sol_i.ye = sol_i.ye(:,keep_indexes);
            
            if ~obs(iobs).IsFinished(sol_i)
                finished = false;
                break
            end
        end
    end
end
