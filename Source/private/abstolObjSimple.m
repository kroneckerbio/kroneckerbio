function abstol = abstolObjSimple(m, con, obj, sol)
% TODO: consider that parameters may be dominated by partial dG/dp
% TODO: consider objective functions that have continuous portions
% TODO: consider multiple objective functions

% Constants
nX = m.nX;
nV = (size(sol.y,1) - nX) / nX; % Will be zero if gradient not computed
nt = length(sol.x);

%% Find the maximum influence point for each parameter

% Initialize the impact on the gradient by each sensitivity
impact = zeros(nX*(1+nV),nt);
dGdx = zeros(nX,nt);

% Cycle through each timepoint and apply dG/dx to get the effect on
for i = 1:nt
    dGdxi = vec(obj.dGdx(sol.x(i), sol, con.u, 1)); % Evaluate dGdx at time ti
    temp = sol.y(:,i); % Pull out the solution at one timepoint
    temp = reshape(temp, nX,1+nV); % Reshape to expose x to dG/dx
    temp = diag(dGdxi) * temp; % Multiply by dG/dx to get the impact on G and dG/dp
    dGdx(:,i) = dGdxi;
    impact(:,i) = temp(:);
end

% Condense to the maximum impact
impact = abs(impact); % The absoluate value is the important value
impact = max(impact, [], 2); % Max over all time
impact = reshape(impact, nX,1+nV); % x_p
impact = max(impact, [], 1); % Max over all species (Final length = 1+p)

%% Find the tolerance of each sensitivity equivalent to the most important reltol

% Initialize the abstol that should be used once multiplied by reltol
abstol = zeros(nX*(1+nV),nt);

% Cycle through each timepoint and calculate the optimum tolerance
for i = 1:nt
    dGdxi = dGdx(:,i); % Fetch dGdx at time ti
    temp = vec(abs(dGdxi.^-1)) * impact; % abstol_ti_xj_pk = impact_xj / dGdx_pk
    abstol(:,i) = temp(:);
end

% Take the conservative estimate (min) and apply for all time
abstol = min(abstol, [], 2); % xp_ keep vectorized

% Prevent zero values from being returned as this will make ode solvers freak
badAbsTol = (abstol <= 0);
if any(badAbsTol)
    abstol(badAbsTol) = min(abstol(~badAbsTol));
end

end