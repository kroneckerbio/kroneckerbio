function h = plotCurvatureExperiment(m, sim, varargin)

if ~isempty(sim.sol.idata)
    % Is a continuous simulation
    d2ydT2 = sim.d2ydT2(sim.t);
else
    % Is a sampled simulation
    d2ydT2 = sim.d2ydT2;
end

h = plot(sim.t, d2ydT2, varargin{:});

xlabel('Time');
ylabel('Curvature');
title(sim.Name);
legend({m.Outputs.Name});
