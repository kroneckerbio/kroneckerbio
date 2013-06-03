function obj = constructObjectiveSS(m)

% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

obj.update      = @update;
obj.g           = @g;
obj.dgdx        = @dgdx;
obj.dgdp        = @dgdp;
obj.G           = @G;
obj.dGdx        = @dGdx;
obj.dGdp        = @dGdp;
obj = orderfields(obj);

scale = 1e1;

%% update function
    function objNew = update(mNew)
        objNew = constructObjectiveSS(mNew);
    end

    function val = g(t,x,u,nSim)
        val = m.f(t, x, u);
        val = scale*val.'*val;
    end
    function val = dgdx(t,x,u,nSim)
        val = scale*2*m.f(t, x, u).'*m.dfdx(t,x,u);
    end
    function val = dgdp(t,x,u,nSim)
        val = scale*2*m.f(t, x, u).'*m.dfdp(t,x,u);
    end
    function [val deltaTimes] = G(t,x,u,nSim)
        val = 0;
        deltaTimes = [];
    end
    function val = dGdx(t,x,u,nSim)
        val = zeros(1, m.nX);
    end
    function val = dGdp(t,x,u,nSim)
        val = zeros(1, m.nP);
    end
end
