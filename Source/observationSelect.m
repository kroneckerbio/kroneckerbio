function obs = observationSelect(tGet)

obs.Type = 'Observation.Select';
obs.Complex = false;

obs.tF = max([0, row(tGet)]);
obs.DiscreteTimes = row(tGet);

obs.Simulation = @simulation;
obs.Sensitivity = @sensitivity;

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        sim.Type = 'Simulation.Select';
        sim.Name = int.Name;
        sim.t    = int.t;
        sim.x    = int.x;
        sim.u    = int.u;
        sim.y    = int.y;
        sim.int  = int;
    end

    function sim = sensitivity(int)
        sim = simulation(int);
        sim.dxdT = int.dxdT;
        sim.dudT = int.dudT;
        sim.dydT = int.dydT;
    end
end
