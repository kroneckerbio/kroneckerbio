function obj = constructObjectiveXgtP(m, xInd, pInd)
% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

% Compute the matrix that selects out the interesting species.
C = eye(m.nX);
C = C(xInd,:);

obj.update = @update;
obj.G = @G;
obj.dGdx = @dGdx;
obj.dGdp = @dGdp;

    function objNew = update(mNew)
        objNew = constructObjectiveXgtP(mNew, xInd, pInd);
    end

    function [val stopTimes] = G(t,xSol,u,nSim)
        nT = length(xSol(nSim).x);
        iM = find(xSol(nSim).y(xInd, :) == max(xSol(nSim).y(xInd, :)), 1, 'first');
        tt = fminbnd(@(t) -deval(xSol(nSim), t, xInd), xSol(nSim).x(max(1, iM-1)), xSol(nSim).x(min(iM+1, nT)));
        tt = tt(1);
        xt = deval(xSol(nSim), tt, 1:m.nX);
        % (2) Computing objective function value
        val       =  m.p(pInd) - C*xt;
        stopTimes =  tt;
    end

    function val = dGdx(t,xSol,u,nSim)
        val = zeros(1, m.nX);
        nT = length(xSol(nSim).x);
        iM = find(xSol(nSim).y(xInd, :) == max(xSol(nSim).y(xInd, :)), 1, 'first');
        tt = fminbnd(@(t) -deval(xSol(nSim), t, xInd), xSol(nSim).x(max(1, iM-1)), xSol(nSim).x(min(iM+1, nT)));
        tt = tt(1);
        if t == tt
            val = -C;
        end
    end

    function val = dGdp(t,xSol,u,nSim)
        val = zeros(1, m.nP);
        nT = length(xSol(nSim).x);
        iM = find(xSol(nSim).y(xInd, :) == max(xSol(nSim).y(xInd, :)), 1, 'first');
        tt = fminbnd(@(t) -deval(xSol(nSim), t, xInd), xSol(nSim).x(max(1, iM-1)), xSol(nSim).x(min(iM+1, nT)));
        tt = tt(1);
        if t == tt
            val(pInd) = 1;
        end
    end
end