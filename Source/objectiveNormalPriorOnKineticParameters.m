function obj = objectiveNormalPriorOnKineticParameters(kbar, Vkbar, name)
%obj = objectiveNormalPriorOnKineticParameters(kbar, Vkbar, name)

% Clean up inputs
if nargin < 3
    name = [];
end

% Check inputs
n = numel(kbar);
assert(ismatrix(Vkbar) && all(size(Vkbar) == [n,n]), 'KroneckerBio:constructObjectiveParameterNormal:Vsize', 'Input "Vkbar" must be a square matrix of numel(kbar)')

if isempty(name)
	name = 'NormalPriorOnKineticParameters';
end

% Objective structure
obj.Type = 'Objective.Information.NormalPriorOnKineticParameters';
obj.Name = name;

obj.Continuous    = false;
obj.Complex       = false;
obj.DiscreteTimes = 0;

obj.G      = @G;
obj.dGdk   = @dGdk;
obj.d2Gdk2 = @d2Gdk2;

obj.p      = @p;
obj.logp   = @logp;
obj.F      = @F;
obj.Fn     = @Fn;

obj = pastestruct(objectiveZero, obj);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [val, stopTimes] = G(int)
        Tk = int.k(int.UseParams);
        Tbark = kbar(int.UseParams);
        VTbark = Vkbar(int.UseParams,int.UseParams);
        FTbark = infoinv(VTbark);
        
        stopTimes = 0;
        diff = Tk - Tbark;
        val = diff.' * FTbark * diff;
    end

    function val = dGdk(int)
        nk = numel(int.k);
        Tk = int.k(int.UseParams);
        Tbark = kbar(int.UseParams);
        VTbark = Vkbar(int.UseParams,int.UseParams);
        FTbark = infoinv(VTbark);
        
        val = sparse(nk,1);
        val(int.UseParams) = 2 * (FTbark * (Tk - Tbark));
    end

    function val = d2Gdk2(int)
        nk = numel(int.k);
        VTbark = Vkbar(int.UseParams,int.UseParams);
        FTbark = infoinv(VTbark);
        
        val = sparse(nk,nk);
        val(int.UseParams,int.UseParams) = 2 * FTbark;
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Information theory %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Likelihood function
    function val = p(sol)
        nk = numel(sol.k);
        Tk = sol.k(sol.UseParams);
        Tbark = kbar(sol.UseParams);
        VTbark = Vkbar(sol.UseParams,sol.UseParams);
        FTbark = infoinv(VTbark);

        diff = Tk - Tbark;
        val = (2*pi).^(-nTk/2) * det(VTbark).^(-1/2) * exp(-1/2 * diff.' * FTbark * diff);
    end

%% Log likelihood
    function val = logp(sol)
        nk = numel(sol.k);
        Tk = sol.k(sol.UseParams);
        Tbark = kbar(sol.UseParams);
        VTbark = Vkbar(sol.UseParams,sol.UseParams);
        FTbark = infoinv(VTbark);

        diff = Tk - Tbark;
        val = (-nTk/2) * log(2*pi) + -1/2 * sum(log(infoeig(VTbark))) + -1/2 * diff.' * FTbark * diff;
    end

%% Fisher information
    function val = F(dxdTSol)
        VTbark = Vkbar(sol.UseParams,sol.UseParams);

        if int.Normalized
            Tbark = kbar(sol.UseParams);
            VTbark = spdiags(Tbark.^(-1),0,nTk,nTk) * VTbark * spdiags(Tbark.^(-1),0,nTk,nTk);
        end
        
        % This objective only provides information on the first nTk parameters
        T = [sol.k(sol.UseParams); sol.s(sol.UseSeeds); sol.q(sol.UseInputControls); sol.h(sol.UseDoseControls)];
        nT = numel(T);
        val = zeros(nT,nT);
        
        val(1:nTk,1:nTk) = infoinv(VTbark);
    end
end
