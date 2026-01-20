function y_ss = iterate_steady_state(der, jac, ic, nx, abstol, reltol, basal_discontinuities, timescale)

% Set up times for first simulation
t0 = 0;
tF = t0 + timescale;

% Ensure we cover the discontinuities
if ~isempty(basal_discontinuities)
    tF = max(tF, max(basal_discontinuities) + timescale);
end

% Repeatedly simulate the timescale interval until we reach steady state
reached_steady_state = false;
while ~reached_steady_state
    
    % Integrate
    sol_tmp = accumulateOdeFwdSimp(der, jac, t0, tF, ic, ...
        basal_discontinuities, t0+timescale, ...
        1:nx, reltol, abstol, [], [], []);
    y_ss = sol_tmp.y(:,end);
    
    % Check for steady state
    abs_diff = abs(y_ss - ic);
    rel_diff = abs_diff ./ abs(ic);
    max_err = max(min(abs_diff - abstol, rel_diff - reltol));
    reached_steady_state = max_err <= 0;
    
    % Set up for next iteration
    ic = y_ss;
    t0 = tF;
    tF = tF + timescale;
    
end

end