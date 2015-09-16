function obj = objectiveWeightedSumOfSquaresNonNeg(outputlist, timelist, sd, measurements, name)
% obj = objectiveWeightedSumOfSquaresNonNeg(outputlist, timelist, sd, measurements, name)
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
assert(numel(timelist) == n, 'KroneckerBio:objectiveWeightedSumOfSquaresNonNeg:timelist', 'Input "timelist" must be a vector length of "outputlist".')
assert(numel(measurements) == n || isempty(measurements), 'KroneckerBio:objectiveWeightedSumOfSquaresNonNeg:measurements', 'Input "measurements" must be a vector length of "outputlist" or empty.')

% Make measurements nonnegative
measurements = max(measurements, 0);
floored = (measurements == 0);

% Find unique timelist
discrete_times = row(unique(timelist));

% Default is real data
if isempty(name)
    name = 'WeightedSumOfSquaresNonNeg';
end

obj.Type = 'Objective.Data.WeightedSumOfSquaresNonNeg';
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

obj = pastestruct(objectiveZero(), obj);

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Helper functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [ybar, sigma] = evaluate_sol(sol)
        nx = sol.nx;

        % Evaluate the ODE solver structure
        if isempty(sol.idata)
            % The solution is already discretized
            ind = lookupmember(discrete_times, sol.x);
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
            ind = lookupmember(discrete_times, sol.x);
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
            ix = x(:,ind);
            idxdT = reshape(dxdT(:,ind), nx,nT); %x_T
            iC1 = sol.C1(outputlist(i),:);
            
            % Compute this y
            y = iC1*ix + sol.C2(outputlist(i),:)*u(:,ind) + sol.c(outputlist(i));
            
            % Compute dy/dT
            dydT(i,:) = iC1 * idxdT;
            
            % Compute expected V matrix
            sigma(i) = sd(timelist(i), outputlist(i), y);
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val, discrete] = G(sol)
        % Evaluate solution
        [ybar, sigma] = evaluate_sol(sol);
        
        % Chi2 is non-floored version
        % Compute weighted square error for each data point
        e = ybar(~floored,1) - measurements(~floored,1);
        squareerror = 2 * log(sigma(~floored,1)) + (e./sigma(~floored,1)).^2;
        
        % Area is floored version
        % Integrate Gaussian from -inf to 0 for each data point
        area = -2 .* log(erfc(ybar(floored,1) ./ (sqrt(2) * sigma(floored,1))));
        
        % Goal function
        val = sum(squareerror) + sum(area);
        
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
            flooredt = floored(indt);
            timelistt = timelist(indt);
            outputlistt = outputlist(indt);
            measurementst = measurements(indt);
            
            % Compute e for each datapoint that has a matching time to t
            ybar = zeros(ylength,1);
            dydx = zeros(ylength,nx);
            sigma = zeros(ylength,1);
            dsigmady = zeros(ylength,1);
            for i = 1:ylength
                dydxi = sol.dydx(timelistt(i), xt, ut);
                dydx(i,:) = dydxi(outputlistt(i),:);
                ybar(i) = yt(outputlistt(i));
                [sigma(i), dsigmady(i)] = sd(timelistt(i), outputlistt(i), ybar(i));
            end
            
            % Non-floored data
            enotfloored = ybar(~flooredt,1) - measurementst(~flooredt,1);
            C1notfloored = dydx(~flooredt,:);
            sigmanotfloored = sigma(~flooredt,1);
            dsigmadynotfloored = dsigmady(~flooredt,1);
            dsigmadxfloored = bsxfun(@times, dsigmadynotfloored, C1notfloored);
            dchi2dx = vec(2 * sum(bsxfun(@times, sigmanotfloored.^-1, dsigmadxfloored),1) + 2 * sum(bsxfun(@times, enotfloored .* sigmanotfloored.^-2, C1notfloored),1) + -2 * sum(bsxfun(@times, enotfloored.^2 .* sigmanotfloored.^-3, dsigmadxfloored),1));
            
            % Floored data
            yfloored = ybar(flooredt,1);
            dydxfloored = dydx(flooredt,:);
            sigmafloored = sigma(flooredt,1);
            dsigmadyfloored = dsigmady(flooredt,1);
            dsigmadxfloored = bsxfun(@times, dsigmadyfloored, dydxfloored);
            
            % The portion that is not a function of x
            constant = 2 .* exp(-yfloored.^2./(2*sigmafloored.^2)) .* sqrt(2/pi) ./ (erfc(yfloored ./ sqrt(2*sigmafloored.^2)) .* sigmafloored.^2); 
            
            % The portion that varies with x
            variable = bsxfun(@times, -yfloored, dsigmadxfloored) + bsxfun(@times, sigmafloored, dydxfloored);
            
            % Combine in appropriate dimensions
            dareadx = bsxfun(@times, constant, variable);
            dareadx = vec(sum(dareadx, 1));
            
            % Total is sum
            val = dchi2dx + dareadx;
        else
            val = zeros(nx, 1);
        end
    end

    function val = d2Gdx2(t, sol)
        nx = size(sol.C1,2);

        error('Error: Not implemented yet.') % Not done yet
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
            val = 2 * C1.' * bsxfun(@times, sigma.^-1, C1) + ...
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
    function p = p(sol)
        % p = prod(tau*sigma^2)^(-1/2) * exp(-1/2 * prod(((ybar-yhat)/sigma)^2)) for non-floored
        %   * prod(1/2 * erfc(ybar / sqrt(2*sigma^2))) for floored
        
        % Evaluate solution
        [ybar, sigma] = evaluate_sol(sol);
        
        % Chi2 is non-floored version
        e = ybar(~floored) - measurements(~floored);
        chi2 = sum((e./sigma(~floored)).^2);
        
        % Area is floored version
        area = 0.5 * erfc(ybar(floored) ./ (sqrt(2) * sigma(floored)));
        
        % Total probability
        p = (2*pi)^(-n/2) * prod(sigma(~floored)).^-1 * exp(-0.5*chi2) * prod(area);
    end

%% Log likelihood
    function val = logp(sol)
        % logp = -n/2 * log(tau) + -sum(log(sigma)) + -1/2 * sum(((ybar-yhat)/sigma)^2) for non-floored
        %      + log(1/2*erfc(ybar/sqrt(2*sigma^2))) for floored
        
        % Evaluate solution
        [ybar, sigma] = evaluate_sol(sol);
        e = ybar - measurements;
        
        % Log likelihood for non-floored and floored data
        val_nonfloored = -(nnz(~floored))/2 * log(2*pi) + -sum(log(sigma(~floored))) + -1/2 * sum((e(~floored) ./ sigma(~floored)).^2);
        val_floored = sum(log(0.5 * erfc(ybar(floored) ./ (sqrt(2) * sigma(floored)))));
        
        val = val_nonfloored + val_floored;
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
    function [objNew, new_measurements] = AddData(sol)
        nx = size(sol.C1,2);
        
        % Evaluate the ODE solver structure
        if isempty(sol.idata)
            % The solution is already discretized
            x = sol.y(1:nx, lookupmember(discrete_times, sol.x));
            u = sol.u(:,lookupmember(discrete_times, sol.x));
        else
            % The complete solution is provided
            x = deval(sol, discrete_times, 1:nx); % x_t
            u = sol.u(discrete_times);
        end
        
        % New measurements
        new_measurements = zeros(n,1);
        for i = 1:n
            t = discrete_times == timelist(i);
            new_measurements(i) = sol.C1(outputlist(i),:) * x(:,t) + sol.C2(outputlist(i),:) * u(:,t) + sol.c(outputlist(i),:);
            new_measurements(i) = new_measurements(i) + randn*sd(timelist(i),outputlist(i),new_measurements(i));
        end
        
        % Non-negative
        new_measurements(new_measurements < 0) = 0;
        
        objNew = objectiveWeightedSumOfSquaresNonNeg(outputlist, timelist, sd, new_measurements, name);
    end
end
