function obj = objectiveLogNormalPriorOnKineticParameters(kbar, Vkbar, name)
%obj = objectiveNormalPriorOnKineticParameters(kbar, Vkbar, name)
% kbar and Vkbar are non-normalized as inputs

% Clean up inputs
if nargin < 3
    name = [];
end

% Check inputs
n = numel(kbar);
assert(ismatrix(Vkbar) && all(size(Vkbar) == [n,n]), 'KroneckerBio:objectiveNormalizedNormalPriorOnKineticParameters:Vsize', 'Input "Vkbar" must be a square matrix of numel(kbar)')

if isempty(name)
	name = 'NormalizedNormalPriorOnKineticParameters';
end

% Objective structure
obj.Type = 'Objective.Information.NormalizedNormalPriorOnKineticParameters';
obj.Name = name;

obj.Continuous    = false;
obj.Complex       = false;
obj.Linked        = 0;
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

    function [val, stopTimes] = G(sol)
        nk = numel(sol.k);
        nTk = nnz(sol.UseParams);
        Tk = sol.k(sol.UseParams);
        Tbark = kbar(sol.UseParams);
        VTbark = Vkbar(sol.UseParams,sol.UseParams);
        
        % Normalize
        logTk = log(Tk);
        logTbark = log(Tbark);
        VlogTbark = spdiags(Tbark.^(-1),0,nTk,nTk) * VTbark * spdiags(Tbark.^(-1),0,nTk,nTk);
        
        FlogTbark = infoinv(VlogTbark);

        stopTimes = 0;
        diff = logTk - logTbark;
        val = diff.' * FlogTbark * diff;
    end

    function val = dGdk(sol)
        nk = numel(sol.k);
        nTk = nnz(sol.UseParams);
        Tk = sol.k(sol.UseParams);
        Tbark = kbar(sol.UseParams);
        VTbark = Vkbar(sol.UseParams,sol.UseParams);
        
        % Normalize
        logTk = log(Tk);
        logTbark = log(Tbark);
        VlogTbark = spdiags(Tbark.^(-1),0,nTk,nTk) * VTbark * spdiags(Tbark.^(-1),0,nTk,nTk);
        
        FlogTbark = infoinv(VlogTbark);

        val = zeros(nk,1);
        val(sol.UseParams) = 2 * (diag(Tk.^(-1)) * (FlogTbark * (logTk - logTbark)));
    end

    function val = d2Gdk2(sol)
        nk = numel(sol.k);
        nTk = nnz(sol.UseParams);
        Tk = sol.k(sol.UseParams);
        Tbark = kbar(sol.UseParams);
        VTbark = Vkbar(sol.UseParams,sol.UseParams);
        
        % Normalize
        logTk = log(Tk);
        logTbark = log(Tbark);
        VlogTbark = spdiags(Tbark.^(-1),0,nTk,nTk) * VTbark * spdiags(Tbark.^(-1),0,nTk,nTk);
        
        FlogTbark = infoinv(VlogTbark);

        val = zeros(nk,nk);
        val(sol.UseParams,sol.UseParams) = 2 * diag(Tk.^(-1)) * FTbark * diag(Tk.^(-1)) - 2 * diag(Tk.^(-2)) * diag(FlogTbark * (logTk - logTbark));
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Information theory %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Likelihood function
    function val = p(sol)
        nk = numel(sol.k);
        nTk = nnz(sol.UseParams);
        Tk = sol.k(sol.UseParams);
        Tbark = kbar(sol.UseParams);
        VTbark = Vkbar(sol.UseParams,sol.UseParams);
        
        % Normalize
        logTk = log(Tk);
        logTbark = log(Tbark);
        VlogTbark = spdiags(Tbark.^(-1),0,nTk,nTk) * VTbark * spdiags(Tbark.^(-1),0,nTk,nTk);
        
        FlogTbark = infoinv(VlogTbark);

        diff = logTk - logTbark;
        val = (2*pi).^(-nTk/2) * det(VlogTbark).^(-1/2) * exp(-1/2 * diff.' * FlogTbark * diff);
    end

%% Log likelihood
    function val = logp(sol)
        nk = numel(sol.k);
        nTk = nnz(sol.UseParams);
        Tk = sol.k(sol.UseParams);
        Tbark = kbar(sol.UseParams);
        VTbark = Vkbar(sol.UseParams,sol.UseParams);
        
        % Normalize
        logTk = log(Tk);
        logTbark = log(Tbark);
        VlogTbark = spdiags(Tbark.^(-1),0,nTk,nTk) * VTbark * spdiags(Tbark.^(-1),0,nTk,nTk);
        
        FlogTbark = infoinv(VlogTbark);

        diff = logTk - logTbark;
        val = (-nTk/2) * log(2*pi) + -1/2 * sum(log(infoeig(VlogTbark))) + -1/2 * diff.' * FlogTbark * diff;
    end

%% Fisher information
    function val = F(sol)
        nTk = nnz(sol.UseParams);
        VTbark = Vkbar(sol.UseParams,sol.UseParams);

        % This objective only provides information on the first nTk parameters
        T = [sol.k(sol.UseParams); sol.s(sol.UseSeeds); sol.q(sol.UseControls)];
        nT = numel(T);
        val = zeros(nT,nT);
        
        val(1:nTk,1:nTk) = infoinv(VTbark);
    end

    function val = Fn(sol)
        nTk = nnz(sol.UseParams);
        Tbark = kbar(sol.UseParams);
        VTbark = Vkbar(sol.UseParams,sol.UseParams);
        
        % Normalize
        VlogTbark = spdiags(Tbark.^(-1),0,nTk,nTk) * VTbark * spdiags(Tbark.^(-1),0,nTk,nTk);

        % This objective only provides information on the first nTk parameters
        T = [sol.k(sol.UseParams); sol.s(sol.UseSeeds); sol.q(sol.UseControls)];
        nT = numel(T);
        val = zeros(nT,nT);
        
        val(1:nTk,1:nTk) = infoinv(VlogTbark);
    end
end
