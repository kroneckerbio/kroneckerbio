function obj = constructObjectiveMultipleOutputsSSE(m, outputIndices, dataArray)

% argument processing and error checking
if nargin<2 || isempty(outputIndices)
    outputIndices = 1;
end

assert(isnumeric(outputIndices), '', 'The second argument to constructObjectiveSpeciesSSE must be an species index.');
assert(outputIndices <= m.nY, '', 'The species index %d is greater than the total number of speciess %d.', outputIndices, m.nY);

% Compute the matrix that selects out the interesting species.
C = eye(m.nY);
C = C(outputIndices,:);

numOutputs = length(dataArray);

for i = 1:numOutputs
    objArray(i) = constructObjectiveOutputSSE(m, outputIndices(i), dataArray(i));
end

obj.G = @G;
obj.dGdx = @dGdx;
obj.dGdp = @dGdp;
obj.d2Gdx2 = @d2Gdx2;
obj.d2Gdp2 = @d2Gdp2;
obj.d2Gdxdp = @d2Gdxdp;
obj.d2Gdpdx = @d2Gdpdx;
obj.update = @update;

%% The objective function
    function objNew = update(mNew)
        objNew = constructObjectiveMultipleOutputsSSE(mNew, outputIndices, dataArray);
    end

    function [val stopTimes] = G(t,xSol,u,nSim)
        val = 0;
        stopTimes = [];
        
        for k = 1:numOutputs
            [vali stopTimesi] = objArray(i).G(t,xSol,u,nSim);
            val = val + vali;
            stopTimes = [stopTimes stopTimesi];
        end
    end

    function val = dGdx(t,xSol,u,nSim)
        val = zeros(1, m.nX);
        for k = 1:numOutputs
            vali = objArray(i).dGdx(t,xSol,u,nSim);
            val = val + vali;
        end
    end

    function val = dGdp(t,xSol,u,nSim)
        val = zeros(1, m.nP);
        for k = 1:numOutputs
            vali = objArray(i).dGdp(t,xSol,u,nSim);
            val = val + vali;
        end
    end

    function val = d2Gdx2(t,xSol,u,nSim)
        val = zeros(m.nX, m.nX);
        for k = 1:numOutputs
            vali = objArray(i).d2Gdx2(t,xSol,u,nSim);
            val = val + vali;
        end
    end

    function val = d2Gdp2(t,xSol,u,nSim)
        val = zeros(m.nP, m.nP);
        for k = 1:numOutputs
            vali = objArray(i).d2Gdp2(t,xSol,u,nSim);
            val = val + vali;
        end
    end

    function val = d2Gdpdx(t,xSol,u,nSim)
        val = zeros(m.nX, m.nP);
        for k = 1:numOutputs
            vali = objArray(i).d2Gdxdp(t,xSol,u,nSim);
            val = val + vali;
        end
    end

    function val = d2Gdxdp(t,xSol,u,nSim)
        val = zeros(m.nP, m.nX);
        for k = 1:numOutputs
            vali = objArray(i).d2Gdpdx(t,xSol,u,nSim);
            val = val + vali;
        end
    end
end