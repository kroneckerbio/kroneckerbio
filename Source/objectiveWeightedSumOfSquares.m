function obj = objectiveWeightedSumOfSquares(outputlist, timelist, sd, measurements, name)
% obj = objectiveWeightedSumOfSquares(outputlist, timelist, sd, measurements, name)
% returns X2, chi-square, the weighted least squares statistic

% sd: @(t,yInd,yVal) returns standard error and gradient given t, index of output, and
% value of output

% Clean up inputs
if nargin < 5
    name = [];
    if nargin < 4
        measurements = [];
    end
end

% Standardize inputs
outputlist = vec(outputlist);
timelist = vec(timelist);
measurements = vec(measurements);

% Check inputs
n = numel(outputlist);
assert(numel(timelist) == n, 'KroneckerBio:objectiveWeightedSumOfSquares:timelist', 'Input "timelist" must be a vector length of "outputlist".')
assert(numel(measurements) == n || isempty(measurements), 'KroneckerBio:objectiveWeightedSumOfSquares:measurements', 'Input "measurements" must be a vector length of "outputlist" or empty.')

% Find unique timelist
discrete_times = row(unique(timelist));

if isempty(name)
    name = 'WeightedSumOfSquares';
end

obj.Type = 'Objective.Data.WeightedSumOfSquares';
obj.Name = name;
obj.DiscreteTimes = discrete_times;
obj.Continuous = false;
obj.Complex = false;

obj.G = @G;
obj.dGdx = @dGdx;
obj.d2Gdx2 = @d2Gdx2;

obj.p = @p;
obj.logp = @logp;
obj.F = @F;
obj.Fn = @Fn;

obj.AddData = @AddData;

obj = pastestruct(objectiveZero, obj);

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Helper functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [ybar, sigma] = evaluate_sol(sol)
        nx = sol.nx;

        % Evaluate the ODE solver structure
        if isempty(sol.idata)
            % The solution is already discretized
            ind = lookup(discrete_times, sol.x);
            x = sol.y(1:nx,ind);
            u = sol.u(:,ind);
            y = sol.y_(discrete_times, x, u);
        else
            % The complete solution is provided
            x = deval(sol, discrete_times, 1:nx); % x_t
            u = sol.u(discrete_times);
            y = sol.y_(discrete_times, x, u);
        end
        
        % Get ybar and sigma for every point evaluated
        ybar = zeros(n,1);
        sigma = zeros(n,1);
        for i = 1:n
            it = (discrete_times == timelist(i));
            ybar(i) = y(outputlist(i), it);
            sigma(i) = sd(timelist(i), outputlist(i), ybar(i));
        end
    end

    function [dydT, sigma] = evaluate_grad(sol)
        T = [sol.k(sol.UseParams); sol.s(sol.UseSeeds); sol.q(sol.UseControls)];
        nT = numel(T);
        nx = size(sol.C1,2);

        % Evaluate the ODE solver structure
        dxdTStart = nx+1;
        dxdTEnd   = nx+nx*nT;
        if isempty(sol.idata)
            % The solution is already discretized
            ind = lookup(discrete_times, sol.x);
            x = sol.y(1:nx,ind);
            u = sol.u(:,ind);
            dxdT = sol.y(dxdTStart:dxdTEnd,ind); % xT_t
        else
            % The complete solution is provided
            joint = deval(sol, discrete_times, 1:dxdTEnd); % x+xT_t
            x = joint(1:nx,:); % x_t
            u = sol.u(discrete_times); % u_t
            dxdT = joint(dxdTStart:dxdTEnd,:); % xT_t
        end
        
        % Get dydT for every point
        dydT = zeros(n, nT);
        sigma = zeros(n,1);
        for i = 1:n
            ind = find(discrete_times == timelist(i), 1);
            it = timelist(i);
            ix = x(:,ind);
            iu = u(:,ind);
            idxdT = reshape(dxdT(:,ind), nx,nT); %x_T
            %iC1 = sol.C1(outputlist(i),:);
            dydx = sol.dydx(timelist(i), ix, iu);
            
            % Compute this y
            y = sol.y_(it, ix, iu);
            %y = iC1*ix + sol.C2(outputlist(i),:)*u(:,ind) + sol.c(outputlist(i));
            
            % Compute dy/dT
            dydT(i,:) = dydx * idxdT;
            
            % Compute expected V matrix
            sigma(i) = sd(it, outputlist(i), y);
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val, discrete] = G(sol)
        % Evaluate solution
        [ybar, sigma] = evaluate_sol(sol);
        e = ybar - measurements;

        % Goal function
        val = sum(2 * log(sigma) + (e./sigma).^2);
        
        % Return discrete times as well
        discrete = discrete_times;
    end

    function val = dGdx(t, sol)
        nx = sol.nx;

        %Find all data points that have a time that matches t
        indt = find(t == timelist);
        ylength = length(indt);
        
        if ylength > 0
            % Evaluate the ODE solver structure
            if isempty(sol.idata)
                % The solution is already discretized
                % Get the solution at the matching timepoint
                tInd = (sol.x == t);
                xt = sol.y(1:nx, tInd);
                ut = sol.u(:,tInd);
                yt = sol.y_(t, xt, ut);
            else
                % The complete solution is provided
                xt = deval(sol, t, 1:nx);
                ut = sol.u(t);
                yt = sol.y_(t, xt, ut);
            end
            
            % Extract the data points with time t
            timelistt = timelist(indt);
            outputlistt = outputlist(indt);
            measurementst = measurements(indt);
            
            % Compute e for each datapoint that has a matching time to t
            ybar = zeros(ylength,1);
            dydx = zeros(ylength,nx); % C1 = dydx
            sigma = zeros(ylength,1);
            dsigmady = zeros(ylength,1);
            for i = 1:ylength
                dydxi = sol.dydx(timelistt(i), xt, ut);
                dydx(i,:) = dydxi(outputlistt(i),:);
                ybar(i) = yt(outputlistt(i));%dydx(i,:) * xt + dydui(outputlistt(i),:) * ut + sol.c(outputlistt(i),:);
                [sigma(i), dsigmady(i)] = sd(timelistt(i), outputlistt(i), ybar(i));
            end
            
            % Gradient function
            e = ybar - measurementst;
            dsigmadx = bsxfun(@times, dsigmady, dydx);
            val = vec(2 * sum(bsxfun(@times, sigma.^-1, dsigmadx),1) + 2 * sum(bsxfun(@times, e .* sigma.^-2, dydx),1) + -2 * sum(bsxfun(@times, e.^2 .* sigma.^-3, dsigmadx),1));
        else
            val = zeros(nx, 1);
        end
    end

    function val = d2Gdx2(t,sol)
        nx = size(sol.C1,2);

        error('Error: Not implemented yet.')
        % THIS STILL NEEDS TO BE VERIFIED
        %Find all data points that have a time that matches t
        ind = find(t == timelist);
        tlength = length(ind);
        if tlength > 0
            C1 = zeros(tlength,nx);
            sigma = zeros(tlength,1);
            dsigmady = zeros(tlength,1);
            d2sigmady2 = zeros(tlength,1);
            for i = 1:tlength
                C1(i,:) = sol.C1(outputlist(ind(i)),:);
                y = C1(i,:)*x + sol.C2(outputlist(ind(i)),:) * u + sol.c(outputlist(ind(i)),:);
                [sigma(i) dsigmady(i) d2sigmady2(i)] = sd(timelist(ind(i)), outputlist(ind(i)), y).^2; %#ok<RHSFN>
            end
            
            % Compute d2Gdx2 = 2*C1'*W*C1 + 4*e'*dWdx*C1 + e'*d2Wdx2*e
            val = 2 * C1.' * bsxfun(@times, sigma.^-1,  C1) + ...
                4 * (bsxfun(@times, e, bsxfun(@times, dsigmady .* sigma.^-3, C1))).' * C1 + ...
                C1.' * bsxfun(@times, e.^2, bsxfun(@times, -2 * (d2sigmady2 .* sigma.^-3 + -3 * dsigmady.^2 .* sigma.^-4), C1));
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
        T = [sol.k(sol.UseParams); sol.s(sol.UseSeeds); sol.q(sol.UseControls)];
        nT = numel(T);

        % Evaluate solution
        [dydT, sigma] = evaluate_grad(sol);
        
        % Assemble variance matrix
        V = spdiags(sigma.^2,0,n,n);

        % Construct normalized sensitivities
        dydT_normalized = dydT * spdiags(T,0,nT,nT); % p along the diagonal

        % Normalized Fisher information matrix
        val = dydT_normalized.' * (V \ dydT_normalized);
        val = symmat(val);
    end

%% AddData
    function objNew = AddData(sol)
        nx = size(sol.C1,2);
        
        % Evaluate the ODE solver structure
        if isempty(sol.idata)
            % The solution is already discretized
            ind = lookup(discrete_times, sol.x);
            x = sol.y(1:nx, ind);
            u = sol.u(:, ind);
        else
            % The complete solution is provided
            x = deval(sol, discrete_times, 1:nx); % x_t
            u = sol.u(discrete_times);
        end
        
        % New measurements
        newMeasurements = zeros(n,1);
        for i = 1:n
            t = discrete_times == timelist(i);
            newMeasurements(i) = sol.C1(outputlist(i),:) * x(:,t) + sol.C2(outputlist(i),:) * u(:,t) + sol.c(outputlist(i),:);
            newMeasurements(i) = newMeasurements(i) + randn * sd(timelist(i), outputlist(i), newMeasurements(i));
        end
        
        objNew = objectiveWeightedSumOfSquares(outputlist, timelist, sd, newMeasurements, name);
    end
end
