function obs = observationEvents(tF, events)
% obs = observationEvents(tF, events)

ne = numel(events);

obs.Complex = false;

obs.tF = tF;
obs.ne = ne;
obs.Events = @event_stitcher;

obs.Simulation = @simulation;
obs.Sensitivity = @sensitivity;
obs.Curvature = @curvature;

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        sim.Type = 'Simulation.System.Events';
        sim.Name = int.Name;
        sim.int  = int;
        sim.ie   = int.ie;
        sim.te   = int.te;
        sim.xe   = int.xe;
        sim.ue   = int.ue;
        sim.ye   = int.ye;
    end

    function sim = sensitivity(int)
        sim = simulation(int);
        sim.Type = 'Simulation.Sensitivity.Events';
        sim.dxedT = int.dxedT;
        sim.duedT = int.duedT;
        sim.dyedT = int.dyedT;
    end

    function sim = curvature(int)
        sim = simulation(int);
        sim.Type = 'Simulation.Curvature.Events';
        sim.d2xedT2 = int.d2xedT2;
        sim.d2uedT2 = int.d2uedT2;
        sim.d2yedT2 = int.d2yedT2;
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
