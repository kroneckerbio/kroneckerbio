function obj = constructObjectiveChiSquareContinuous(m, outputs, points, sd)

% Gut m
temp = m;
clear m
m.nX = temp.nX;
m.nY = temp.nY;
m.nP = temp.nP;
m.c = temp.c;
clear temp

% Constants
nx = m.nX;
ny = m.nY;
nAy = numel(outputs); % Number of active outputs

obj.update = @update;
obj.G = [];
obj.dGdx = [];
obj.dGdp = [];
obj.d2Gdx2 = [];
obj.d2Gdp2 = [];
obj.d2Gdpdx = [];
obj.d2Gdxdp = [];
obj.n = points;
obj.H = [];
obj.Hn = [];
obj.EH = [];
obj.EHn = [];
obj.F = [];
obj.Fn = @Fn;
%obj.FnSimbio = @FnSimbio;
obj.EG = [];
obj.AddData = [];
obj.AddSimulationData = [];
obj.ModelPValue = [];

    function objNew = update(mNew)
        objNew = constructObjectiveChiSquareContinuous(mNew, outputs, points, sd);
    end

    function val = Fn(dxdTSol, T)
        nT = numel(T);
        
        % Evaluate the ODE solver structure
        t     = dxdTSol.x; % _t
        nt    = numel(t);
        C     = dxdTSol.C; % y_x
        nx    = size(C,2);
        x     = dxdTSol.y(1:nx,:); % x_t
        dxdT  = dxdTSol.y((nx+1):(nx*(nT+1)),:); % xT_t
        
        % Convert species to outputs
        y = C * x; % y_t
        dxdT = reshape(dxdT, nx,nT*nt); % x_Tt
        dydT = C * dxdT; % y_Tt
        dydT = reshape(dydT, ny,nT,nt); % y_T_t
        
        % Compute weight
        w = zeros(nAy, nt); % y_t
        for i = 1:nAy
            w(i,:) = sd(t, outputs(i), y(outputs(i),:)).^-2;
        end
        
        % Construct normalization matrix
        dTdlnT = spdiags(T,0,sparse(nT,nT)); % T along the diagonal

        % Create integrand
        integrand = zeros(nT,nT,nt); % T_T_t
        for i = 1:nt
            integrand(:,:,i) = dydT(outputs,:,i).' * spdiags(w(:,i),0,sparse(nAy,nAy)) * dydT(outputs,:,i);
        end
        
        % Integrate by trapezoidal rule
        val = trapz(t, integrand, 3);
        
        % Average information of a single point
        val = val / (ny * (t(end) - t(1)));
        
        % Multiply by the number of points
        val = val * points;
        
        % Symmetrize
        val = symmat(dTdlnT * val * dTdlnT);
    end

%     function val = FnSimbio(simData, p)
%         nP = (size(simData.Data, 2) - m.nX) / m.nX;
%         
%         % Extract x and dxdp
%         times = simData.Time;
%         nt = length(times);
%         x = simData.Data(:,1:nX).'; % x_t
%         dxdp = simData.Data(:,(nX+1):end).'; % xp_t
%         
%         % Convert species to outputs
%         y = C*x; % y_t
%         dxdp = reshape(dxdp, nX,nP*nt); % x_pt
%         dydp = C*dxdp; % y_pt
%         dydp = reshape(dydp, nY,nP,nt); % y_p_t
%         
%         % Compute weight
%         w = zeros(nY, nt); % y_t
%         for i = 1:nY
%             w(i,:) = sd(times, outputs(i), y(i,:)).^-2 * points;
%         end
%         
%         % Construct normalization matrix
%         dpdlnp = spdiags(p,0,sparse(nP,nP)); % p along the diagonal
% 
%         % Create integrand
%         integrand = zeros(nP, nP, nt); % p_p_t
%         for i = 1:nt
%             dydlnp = dydp(:,:,i) * dpdlnp; % y_p
%             integrand(:,:,i) = dydlnp.' * spdiags(w(:,i),0,sparse(nY,nY)) * dydlnp;
%         end
%         
%         % Integrate by trapezoidal rule
%         val = trapz(times, integrand, 3);
%         
%         % Divide by number of species and length of time
%         val = val / (nY * (times(end) - times(1))) ;
%         
%         % Symmetrize
%         val = symmat(val);
%     end
end