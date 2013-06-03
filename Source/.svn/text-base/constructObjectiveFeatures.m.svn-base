function obj = constructObjectiveFeatures(m, outputIndex, data, tOn)
% this function computes the SSE of amplitudes AND timings
% 
% data.tM
% data.tm
% data.yM
% data.ym

nO = length(outputIndex);
for i = 1:nO
    % get indices of the appropriate events
    ii = find(data.ie == i);
    
%     ye = data.ye(i,ii);
%     if ye(2) > ye(1)
%         ii = ii(2:end);
%     end
    
    % get data structures
    data_t.i  = data.ie(ii);
    data_t.t  = data.te(ii);
    data_t.order = data.order;
    data_y.i  = data.ie(ii);
    data_y.t  = data.te(ii);
    data_y.y  = data.ye(i,ii);
    data_y0.t = tOn;
    data_y0.y = data.y0;
    data_y.order = data.order;
    
    % get objective functions
    objv{(i-1)*2 + 1} = constructObjectiveTimingOutputSSE(m, outputIndex(i), data_t, i, tOn);
    objv{(i-1)*2 + 2} = constructObjectiveAmplitudeOutputSSE(m, outputIndex(i), data_y, i, tOn);
end
objv{end+1} = constructObjectiveOutputSSE(m, outputIndex, data_y0);

obj.update      = @update;
obj.g           = @g;
obj.dgdx        = @dgdx;
obj.dgdp        = @dgdp;
obj.G           = @G;
obj.dGdx        = @dGdx;
obj.dGdp        = @dGdp;
obj = orderfields(obj);

%% update function
    function objNew = update(mNew)
        objNew = constructObjectiveFeatures(mNew, outputIndex, data, tOn);
    end

    function [val stopTimes] = G(t, xSol, u, nSim)
        val = 0;
        stopTimes = [];
        for i = 1:length(objv)
            [tmp_val(i) tmp_stopTimes] = objv{i}.G(t, xSol, u, nSim);
            stopTimes = unique([stopTimes tmp_stopTimes]);
        end
        val = sum(tmp_val);
    end

    function val = dGdx(t, xSol, u, nSim)
        val = zeros(1, m.nX);
        for i = 1:length(objv)
            tmp_val = objv{i}.dGdx(t, xSol, u, nSim);
            val = val + tmp_val;
        end
    end

    function val = dGdp(t, xSol, u, nSim)
        val = zeros(1, m.nP);
        for i = 1:length(objv)
            tmp_val = objv{i}.dGdp(t, xSol, u, nSim);
            val = val + tmp_val;
        end
    end

    function val = g(t, x, u, nSim)
        val = 0;
    end

    function val = dgdx(t, x, u, nSim)
        val = zeros(1, m.nX);
    end

    function val = dgdp(t, x, u, nSim)
        val = zeros(1, m.nP);
    end
end