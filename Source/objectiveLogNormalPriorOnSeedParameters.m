function obj = objectiveLogNormalPriorOnSeedParameters(sbar, Vsbar, name)
%obj = objectiveLogNormalPriorOnSeedParameters(m, sbar, Vsbar, name)
% sbar and Vsbar are non-normalized as inputs

% Clean up inputs
if nargin < 3
    name = [];
end

% Check inputs
n = numel(sbar);
assert(ismatrix(Vsbar) && all(size(Vsbar) == [n,n]), 'KroneckerBio:objectiveNormalizedNormalPriorOnSeedParameters:Vsize', 'Input "Vsbar" must be a square matrix of numel(sbar)')

if isempty(name)
	name = 'NormalizedNormalPriorOnSeedParameters';
end

% Objective structure
obj.Type = 'Objective.Information.NormalizedNormalPriorOnSeedParameters';
obj.Name = name;

obj.Continuous    = false;
obj.Complex       = false;
obj.Linked        = 0;
obj.DiscreteTimes = 0;

obj.G      = @G;
obj.dGds   = @dGds;
obj.d2Gds2 = @d2Gds2;

obj.p      = @p;
obj.logp   = @logp;
obj.F      = @F;
obj.Fn     = @Fn;

obj = pastestruct(objectiveZero, obj);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [val stopTimes] = G(int)
        ns = numel(int.s);
        nTs = nnz(int.UseSeeds);
        Ts = int.s(int.UseSeeds);
        Tbars = sbar(int.UseSeeds);
        VTbars = Vsbar(int.UseSeeds,int.UseSeeds);
        
        % Normalize
        logTs = log(Ts);
        logTbars = log(Tbars);
        VlogTbars = spdiags(Tbars.^(-1),0,nTs,nTs) * VTbars * spdiags(Tbars.^(-1),0,nTs,nTs);
        
        FlogTbars = infoinv(VlogTbars);

        stopTimes = 0;
        diff = logTs - logTbars;
        val = diff.' * FlogTbars * diff;
    end

    function val = dGds(int)
        ns = numel(int.s);
        nTs = nnz(int.UseSeeds);
        Ts = int.s(int.UseSeeds);
        Tbars = sbar(int.UseSeeds);
        VTbars = Vsbar(int.UseSeeds,int.UseSeeds);
        
        % Normalize
        logTs = log(Ts);
        logTbars = log(Tbars);
        VlogTbars = spdiags(Tbars.^(-1),0,nTs,nTs) * VTbars * spdiags(Tbars.^(-1),0,nTs,nTs);
        
        FlogTbars = infoinv(VlogTbars);

        val = sparse(ns,1);
        val(int.UseSeeds) = 2 * (diag(Ts.^(-1)) * (FlogTbars * (logTs - logTbars)));
    end

    function val = d2Gds2(int)
        ns = numel(int.s);
        nTs = nnz(int.UseSeeds);
        Ts = int.s(int.UseSeeds);
        Tbars = sbar(int.UseSeeds);
        VTbars = Vsbar(int.UseSeeds,int.UseSeeds);
        
        % Normalize
        logTs = log(Ts);
        logTbars = log(Tbars);
        VlogTbars = spdiags(Tbars.^(-1),0,nTs,nTs) * VTbars * spdiags(Tbars.^(-1),0,nTs,nTs);
        
        FlogTbars = infoinv(VlogTbars);

        val = sparse(ns,ns);
        val(int.UseSeeds,int.UseSeeds) = 2 * diag(Ts.^(-1)) * FlogTbars * diag(Ts.^(-1)) - 2 * diag(Ts.^(-2)) * diag(FlogTbars * (logTs - logTbars));
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Information theory %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Likelihood function
    function val = p(sol)
        ns = numel(sol.s);
        nTs = nnz(sol.UseSeeds);
        Ts = sol.s(sol.UseSeeds);
        Tbars = sbar(sol.UseSeeds);
        VTbars = Vsbar(sol.UseSeeds,sol.UseSeeds);
        
        % Normalize
        logTs = log(Ts);
        logTbars = log(Tbars);
        VlogTbars = spdiags(Tbars.^(-1),0,nTs,nTs) * VTbars * spdiags(Tbars.^(-1),0,nTs,nTs);
        
        FlogTbars = infoinv(VlogTbars);

        diff = logTs - logTbars;
        val = (2*pi).^(-nTs/2) * det(VlogTbars).^(-1/2) * exp(-1/2 * diff.' * FlogTbars * diff);
    end

%% Log likelihood
    function val = logp(sol)
        ns = numel(sol.s);
        nTs = nnz(sol.UseSeeds);
        Ts = sol.s(sol.UseSeeds);
        Tbars = sbar(sol.UseSeeds);
        VTbars = Vsbar(sol.UseSeeds,sol.UseSeeds);
        
        % Normalize
        logTs = log(Ts);
        logTbars = log(Tbars);
        VlogTbars = spdiags(Tbars.^(-1),0,nTs,nTs) * VTbars * spdiags(Tbars.^(-1),0,nTs,nTs);
        
        FlogTbars = infoinv(VlogTbars);

        diff = logTs - logTbars;
        val = (-nTs/2) * log(2*pi) + -1/2 * sum(log(infoeig(VlogTbars))) + -1/2 * diff.' * FlogTbars * diff;
    end

%% Fisher information
    function val = F(sol)
        nTk = nnz(sol.UseParams);
        nTs = nnz(sol.UseSeeds);
        VTbars = Vsbar(sol.UseSeeds,sol.UseSeeds);

        if sol.Normalized
            Tbars = sbar(sol.UseSeeds);
            VTbars = spdiags(Tbars.^(-1),0,nTs,nTs) * VTbars * spdiags(Tbars.^(-1),0,nTs,nTs);
        end
        
        % This objective only provides information on the first nTs parameters
        T = [sol.k(sol.UseParams); sol.s(sol.UseSeeds); sol.q(sol.UseInputControls); sol.h(sol.UseDoseControls)];
        nT = numel(T);
        val = zeros(nT,nT);
        
        inds = nTk+1:nTk+nTs;
        val(inds,inds) = infoinv(VTbars);
    end
end
