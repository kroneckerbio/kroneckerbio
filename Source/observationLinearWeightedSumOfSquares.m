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
obs.DiscreteTimes = row(unique(timelist));

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
    end

    function obj = objective(measurements)
        assert(numel(measurements) == n , 'KroneckerBio:observationWeightedSumOfSquares:measurements', 'Input "measurements" must be a vector length of "outputlist"')
        obj = objectiveLinearWeightedSumOfSquares(outputlist, timelist, measurements, sd, name);
    end
end

function obj = objectiveLinearWeightedSumOfSquares(outputlist, timelist, measurements, sd, name)
% Find unique timelist
n = numel(outputlist);
discrete_times = row(unique(timelist));

obj.Type = 'Objective.Data.WeightedSumOfSquares';
obj.Name = name;
obj.DiscreteTimes = discrete_times;
obj.Continuous = false;
obj.Complex = false;

obj.tF = max(discrete_times);

obj.G = @G;
obj.dGdx = @dGdx;
obj.d2Gdx2 = @d2Gdx2;

obj.p = @p;
obj.logp = @logp;
obj.F = @F;
obj.Fn = @Fn;

obj.AddData = @AddData;

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

    function val = dGdx(t, int)
        nx = int.nx;
        ny = int.ny;
        
        % Find all data points that have a time that matches t
        ind_t = find(t == timelist);
        n_current = numel(ind_t);
        
        if n_current > 0
            % Extract integration for this time point
            yt = int.y(:,int.t == t);
            dydxt = reshape(int.dydx(:,int.t == t), ny,nx);
            
            % Extract the data points with time t
            timelist_t = timelist(ind_t);
            outputlist_t = outputlist(ind_t);
            measurements_t = measurements(ind_t);
            
            % Compute e for each datapoint that has a matching time to t
            ybar_t = zeros(n_current,1);
            dydx_t = zeros(n_current,nx);
            sigma_t = zeros(n_current,1);
            dsigmady_t = zeros(n_current,1);
            for i = 1:n_current
                ybar_t(i) = yt(outputlist_t(i));
                dydx_t(i,:) = dydxt(outputlist_t(i),:);
                [sigma_t(i), dsigmady_t(i)] = sd(timelist_t(i), outputlist_t(i), ybar_t(i));
            end
            
            % Gradient function
            e = ybar_t - measurements_t;
            dsigmadx = bsxfun(@times, dsigmady_t, dydx_t);
            val = vec(2 * sum(bsxfun(@times, sigma_t.^-1, dsigmadx),1) + 2 * sum(bsxfun(@times, e .* sigma_t.^-2, dydx_t),1) + -2 * sum(bsxfun(@times, e.^2 .* sigma_t.^-3, dsigmadx),1));
        else
            val = zeros(nx, 1);
        end
    end

    function val = d2Gdx2(t, int)
        nx = int.nx;
        
        error()
        % THIS STILL NEEDS TO BE VERIFIED
        % Find all data points that have a time that matches t
        ind_t = find(t == timelist);
        n_current = numel(ind_t);
        
        if n_current > 0
            % Extract integration for this time point
            yt = int.y(:,int.t == t);
            dydxt = reshape(int.y(:,int.t == t), ny,nx);

            ybar_t = zeros(n_current,1);
            dydx_t = zeros(n_current,nx);
            sigma_t = zeros(n_current,1);
            dsigmady_t = zeros(n_current,1);
            d2sigmady2_t = zeros(n_current,1);
            for i = 1:n_current
                ybar_t(i) = yt(outputlist_t(i));
                dydx_t(i,:) = dydxt(outputlist_t(i),:);
                [sigma_t(i), dsigmady_t(i), d2sigmady2_t(i)] = sd(timelist(ind_t(i)), outputlist(ind_t(i)), ybar_t(i));
            end
            
            % Compute d2Gdx2 = 2*C1'*W*C1 + 4*e'*dWdx*C1 + e'*d2Wdx2*e
            val = 2 * dydx_t.' * bsxfun(@times, sigma_t.^-1,  dydx_t) + ...
                4 * (bsxfun(@times, e, bsxfun(@times, dsigmady_t .* sigma_t.^-3, dydx_t))).' * dydx_t + ...
                dydx_t.' * bsxfun(@times, e.^2, bsxfun(@times, -2 * (d2sigmady2_t .* sigma_t.^-3 + -3 * dsigmady_t.^2 .* sigma_t.^-4), dydx_t));
        else
            val = zeros(nx, nx);
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
