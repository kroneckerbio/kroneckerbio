function obs = observationEvents(tF, events)

ne = numel(events);

obs.Complex = false;

obs.tF = tF;
obs.ne = ne;
obs.Events = @event_stitcher;

obs.Simulation = @simulation;
obs.Sensitivity = @sensitivity;

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        sim.Type = 'Simulation.Events';
        sim.Name = int.Name;
        sim.ie   = int.ie;
        sim.te   = int.te;
        sim.xe   = int.xe;
        sim.ue   = int.ue;
        sim.ye   = int.ye;
        sim.int  = int;
    end

    function sim = sensitivity(int)
        sim = simulation(int);
        sim.dxedT = int.dxedT;
        sim.duedT = int.duedT;
        sim.dyedT = int.dyedT;
    end

    function [value, isterminal, direction] = event_stitcher(t,y)
        value = zeros(ne,1);
        isterminal = false(ne,1);
        direction = zeros(ne,1);
        for ie = 1:ne
            value(ie) = events(ie).e(t,y);
            direction(ie) = events(ie).direction;
        end
    end
end
