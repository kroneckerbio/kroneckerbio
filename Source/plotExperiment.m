function h = plotExperiment(m, sim, varargin)

if ~isempty(sim.sol.idata)
    % Is a continuous simulation
    y = sim.y(sim.t);
else
    % Is a sampled simulation
    y = sim.y;
end

h = plot(sim.t, y, varargin{:});

xlabel('Time');
ylabel('Amount');
title(sim.Name);
legend({m.Outputs.Name});
