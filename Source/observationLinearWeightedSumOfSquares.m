function obs = observationLinearWeightedSumOfSquares(outputlist, timelist, sd, name)
% Clean up inputs
if nargin < 4
    name = '';
end

% Standardize inputs
outputlist = vec(outputlist);
timelist = vec(timelist);

% Find unique timelist
discrete_times = row(unique(timelist));

if isempty(name)
    name = 'WeightedSumOfSquares';
end

n = numel(outputlist);
assert(numel(timelist) == n, 'KroneckerBio:objectiveWeightedSumOfSquares:timelist', 'Input "timelist" must be a vector length of "outputlist".')

obs.Type = 'Observation.Data.LinearWeightedSumOfSquares';
obs.Name = name;
obs.Complex = false;

obs.tF = max(timelist);
obs.DiscreteTimes = discrete_times;

obs.Simulation = @simulation;
obs.Objective  = @objective;

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        y_all = int.y;
        
        ybar = zeros(n,1);
        for i = 1:n
            ind = discrete_times == timelist(i);
            ybar(i) = y_all(outputlist(i),ind);
            ybar(i) = ybar(i) + randn * sd(timelist(i), outputlist(i), ybar(i));
        end
        
        sim.Type = 'Simulation.Data.LinearWeightedSumOfSquares';
        sim.Name = name;
        sim.t = timelist;
        sim.y = ybar;
        sim.more.outputlist = outputlist;
        sim.more.sd = sd;
        sim.int = int;
    end

    function obj = objective(measurements)
        measurements = vec(measurements);
        assert(numel(measurements) == n , 'KroneckerBio:observationWeightedSumOfSquares:measurements', 'Input "measurements" must be a vector length of "outputlist"')
        obj = objectiveLinearWeightedSumOfSquares(outputlist, timelist, measurements, sd, name);
    end
end

function obj = objectiveLinearWeightedSumOfSquares(outputlist, timelist, measurements, sd, name)
% Find unique timelist
n = numel(outputlist);
discrete_times = row(unique(timelist));

% Inherit observation
obj = observationLinearWeightedSumOfSquares(outputlist, timelist, sd, name);

obj.Type = 'Objective.Data.WeightedSumOfSquares';
obj.Continuous = false;

obj.G = @G;
obj.dGdy = @dGdy;
obj.d2Gdy2 = @d2Gdy2;

obj.p = @p;
obj.logp = @logp;
obj.F = @F;
obj.Fn = @Fn;

obj = pastestruct(objectiveZero(), obj);

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Helper functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [ybar, sigma] = evaluate_sol(int)
        y_all = int.y;
        
        ybar = zeros(n,1);
        sigma = zeros(n,1);
        for i = 1:n
            ind = discrete_times == timelist(i);
            ybar(i) = y_all(outputlist(i),ind);
            sigma(i) = sd(timelist(i), outputlist(i), ybar(i));
        end
    end

    function [dydT, sigma] = evaluate_grad(int)
        nT = int.nT;
        nx = int.nx;
        
        x_all = int.x;
        u_all = int.u;
        y_all = int.y;
        dxdT_all = int.dxdT;

        % Get dydT for every point
        dydT = zeros(n, nT);
        sigma = zeros(n,1);
        for i = 1:n
            ind = find(discrete_times == timelist(i), 1);
            t_i = timelist(i);
            x_i = x_all(:,ind);
            u_i = u_all(:,ind);
            dxdT_i = reshape(dxdT_all(:,ind), nx,nT); %x_T
            dydx_i = int.dydx(timelist(i), x_i, u_i);
            dydx_i = dydx_i(outputlist(i),:);
            
            % Compute this y
            y_i = y_all(outputlist(i),ind);
            
            % Compute dy/dT
            dydT(i,:) = dydx_i * dxdT_i;
            
            % Compute expected V matrix
            sigma(i) = sd(t_i, outputlist(i), y_i);
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val, discrete] = G(int)
        % Evaluate solution
        [ybar, sigma] = evaluate_sol(int);
        e = ybar - measurements;
        
        % Goal function
        val = sum(2 * log(sigma) + (e./sigma).^2);
        
        % Return discrete times as well
        discrete = discrete_times;
    end

    function val = dGdy(t, int)
        % dGdy = 2 * sigma^-1 * dsigmady + 2 * (sigma - (y-yhat)*dsigmady) * sigma^-2
        ny = int.ny;
        
        % Find all data points that have a time that matches t
        ind_t = find(t == timelist);
        n_current = numel(ind_t);
        
        if n_current > 0
            % Extract integration for this time point
            yt = int.y(:,int.t == t);
            
            % Extract the data points with time t
            timelist_t = timelist(ind_t);
            outputlist_t = outputlist(ind_t);
            measurements_t = measurements(ind_t);
            
            % Compute e for each datapoint that has a matching time to t
            ybar_t = zeros(n_current,1);
            sigma_t = zeros(n_current,1);
            dsigmady_t = zeros(n_current,1);
            for i = 1:n_current
                ybar_t(i) = yt(outputlist_t(i));
                [sigma_t(i), dsigmady_t(i)] = sd(timelist_t(i), outputlist_t(i), ybar_t(i));
            end
            
            % Gradient value
            e = ybar_t - measurements_t; % Y_
            dGdybar = 2 ./ sigma_t .* dsigmady_t + 2 .* e .* (sigma_t - e.*dsigmady_t) ./ sigma_t.^3; % Y_
            val = accumarray(outputlist_t, dGdybar, [ny,1]); % sum the entries associated with the same output
        else
            val = zeros(ny, 1);
        end
    end

    function val = d2Gdy2(t, int)
        % d2Gdy2 = -2 * sigma^-2 * dsigmady + 2 * sigma^-1 * d2sigmady2 + 2 * ((sigma - (y-yhat)*dsigmady) * sigma^-3 * (1 - 3*(y-yhat)*sigma^-1) - (y-yhat)^2 * d2sigmady2)
        ny = int.ny;

        % Find all data points that have a time that matches t
        ind_t = find(t == timelist);
        n_current = numel(ind_t);
        
        if n_current > 0
            % Extract integration for this time point
            yt = int.y(:,int.t == t);

            % Extract the data points with time t
            timelist_t = timelist(ind_t);
            outputlist_t = outputlist(ind_t);
            measurements_t = measurements(ind_t);
            
            ybar_t = zeros(n_current,1);
            sigma_t = zeros(n_current,1);
            dsigmady_t = zeros(n_current,1);
            d2sigmady2_t = zeros(n_current,1);
            for i = 1:n_current
                ybar_t(i) = yt(outputlist_t(i));
                [sigma_t(i), dsigmady_t(i), d2sigmady2_t(i)] = sd(timelist_t(i), outputlist_t(i), ybar_t(i));
            end
            
            % Curvature value
            e = ybar_t - measurements_t;
            d2Gdybar2 = -2 ./ sigma_t.^2 .* dsigmady_t + 2 ./ sigma_t .* d2sigmady2_t + ...
                2 .* ((sigma_t - e .* dsigmady_t) ./ sigma_t.^3 .* (1 - 3 .* e ./ sigma_t) - e.^2 .* d2sigmady2_t);
            val = accumarray(outputlist_t, d2Gdybar2, [ny,1]); % y_ % sum the entries associated with the same output
            val = spdiags(val, 0, ny, ny); % y_y
        else
            val = zeros(ny, ny);
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Information theory %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probability density function
    function val = p(sol)
        % p = tau^(-n/2) * prod(sigma)^-1 * exp(-1/2 * sum(((ybar-yhat)/sigma)^2)
        
        % Evaluate solution
        [ybar, sigma] = evaluate_sol(sol);
        e = ybar - measurements;
        
        val = (2*pi)^(-n/2) * prod(sigma).^-1 * exp(-1/2 * sum((e ./ sigma).^2));
    end

%% Log likelihood
    function val = logp(sol)
        % logp = -n/2 * log(tau) + -sum(log(sigma)) + -1/2 * sum(((ybar-yhat)/sigma)^2)

        % Evaluate solution
        [ybar, sigma] = evaluate_sol(sol);
        e = ybar - measurements;
        
        val = -n/2 * log(2*pi) + -sum(log(sigma)) + -1/2 * sum((e ./ sigma).^2);
    end

%% Fisher information matrix
    function val = F(sol)
        % Evaluate solution
        [dydT, sigma] = evaluate_grad(sol);
        
        % Assemble variance matrix
        V = spdiags(sigma.^2,0,n,n);
        
        % Fisher information matrix
        val = dydT.' * (V \ dydT);
        val = symmat(val);
    end

%% Normalized Fisher information matrix
    function val = Fn(sol)
        T = [sol.k(sol.UseParams); sol.s(sol.UseSeeds); sol.q(sol.UseInputControls); sol.h(sol.UseDoseControls)];
        nT = numel(T);

        % Evaluate solution
        [dydT, sigma] = evaluate_grad(sol);
        
        % Assemble variance matrix
        V = spdiags(sigma.^2,0,n,n);

        % Construct normalized sensitivities
        dydT_normalized = dydT * spdiags(T,0,nT,nT); % T along the diagonal

        % Normalized Fisher information matrix
        val = dydT_normalized.' * (V \ dydT_normalized);
        val = symmat(val);
    end
end
