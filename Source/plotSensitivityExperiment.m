function h = plotSensitivityExperiment(m, sim, varargin)

if ~isempty(sim.sol.idata)
    % Is a continuous simulation
    dydT = sim.dydT(sim.t);
else
    % Is a sampled simulation
    dydT = sim.dydT;
end

h = plot(sim.t, dydT, varargin{:});

xlabel('Time');
ylabel('Sensitivity');
title(sim.Name);
legend({m.Outputs.Name});
