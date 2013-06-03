function obj = constructObjectiveMultipleSpeciesSSE(m, speciesIndices, dataArray)

% argument processing and error checking
if nargin<2 || isempty(speciesIndices)
    speciesIndices = 1;
end

assert(isnumeric(speciesIndices), '', 'The second argument to constructObjectiveSpeciesSSE must be an species index.');
assert(size(data.x, 1) == length(speciesIndices), '', 'The length of speciesIndices must be the same as the number of species in data.x.');
assert(speciesIndices <= m.nX, '', 'The species index %d is greater than the total number of speciess %d.', speciesIndices, m.nX);

% Compute the matrix that selects out the interesting species.
C = eye(m.nX);
C = C(speciesIndices,:);

numSpecies = length(dataArray);

for i = 1:numSpecies
    objArray(i) = constructObjectiveSpeciesSSE(m, speciesIndices(i), dataArray(i));
end

obj.G = @G;
obj.dGdx = @dGdx;
obj.dGdp = @dGdp;

%% The objective function
    function objNew = update(mNew)
        objNew = constructObjectiveMultipleSpeciesSSE(mNew, speciesIndices, dataArray);
    end

    function [val stopTimes] = G(t,xSol,u,nSim)
        val = 0;
        stopTimes = [];
        
        for i = 1:numSpecies
            [vali stopTimesi] = objArray(i).G(t,xSol,u,nSim);
            val = val + vali;
            stopTimes = [stopTimes stopTimesi];
        end
    end

    function val = dGdx(t,xSol,u,nSim)
        for i = 1:numSpecies
            vali = objArray(i).dGdx(t,xSol,u,nSim);
            val = val + vali;
        end
    end

    function val = dGdp(t,xSol,u,nSim)
        for i = 1:numSpecies
            vali = objArray(i).dGdp(t,xSol,u,nSim);
            val = val + vali;
        end
    end
end