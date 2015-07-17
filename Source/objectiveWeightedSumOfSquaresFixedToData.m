function obj = objectiveWeightedSumOfSquaresFixedToData(outputlist, timelist, sd, measurements, name)
% obj = objectiveWeightedSumOfSquaresFixedToData(m, outputlist, timelist, sd, measurements, name)
% returns X2, chi-square, the weighted least squares statistic

% sd: @(t,yInd,yVal) returns standard error given t, index of output, and
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
assert(numel(timelist) == n, 'KroneckerBio:objectiveWeightedSumOfSquaresFixedToData:timelist', 'Input "timelist" must be a vector length of "outputlist".')
assert(numel(measurements) == n || isempty(measurements), 'KroneckerBio:objectiveWeightedSumOfSquaresFixedToData:measurements', 'Input "measurements" must be a vector length of "outputlist" or empty.')

% Find unique timelist
discreteTimes = row(unique(timelist));

% Compute covariance matrix
if ~isempty(measurements)
    sigma = zeros(n,1);
    for j = 1:n
        sigma(j) = sd(timelist(j), outputlist(j), measurements(j));
    end
    V = spdiags(sigma.^2, 0, n, n);
end

if isempty(name)
    name = 'WeightedSumOfSquaresFixedToData';
end

obj.Type = 'Objective.Data.WeightedSumOfSquaresFixedToData';
obj.Name = name;
obj.DiscreteTimes = discreteTimes;
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [val discrete] = G(sol)
        nx = size(sol.C1,2);

        % Evaluate the ODE solver structure
        if isempty(sol.idata)
            % The solution is already discretized
            ind = lookup(discreteTimes, sol.x);
            x = sol.y(1:nx,ind);
            u = sol.u(:,ind);
        else
            % The complete solution is provided
            x = deval(sol, discreteTimes, 1:nx); % x_t
            u = sol.u(discreteTimes);
        end
        
        % Get e for every point evaluated
        e = zeros(n,1);
        for i = 1:n
            t = (discreteTimes == timelist(i));
            e(i,1) = sol.C1(outputlist(i),:) * x(:,t) + sol.C2(outputlist(i),:) * u(:,t) + sol.c(outputlist(i),:) - measurements(i);
        end
        
        % Compute chi-square
        val = e'*(V\e);
        
        if nargout >= 2
            discrete = discreteTimes;
        end
    end

    function val = dGdx(t,sol)
        nx = size(sol.C1,2);

        % Find all data points that have a time that matches t
        ind = find(t == timelist);
        ylength = length(ind);
        if ylength > 0
            % Evaluate the ODE solver structure
            if isempty(sol.idata)
                % The solution is already discretized
                % Get the solution at the matching timepoint
                tInd = (sol.x == t);
                x = sol.y(1:nx, tInd);
                u = sol.u(:,tInd);
            else
                % The complete solution is provided
                x = deval(sol, t, 1:nx);
                u = sol.u(t);
            end
            
            % Compute e for each datapoint that has a matching time to t
            e = zeros(ylength,1);
            C1 = zeros(ylength,nx);
            for i = 1:ylength
                C1(i,:) = sol.C1(outputlist(ind(i)),:);
                e(i,1) = C1(i,:)*x + sol.C2(outputlist(ind(i)),:) * u + sol.c(outputlist(ind(i)),:) - measurements(ind(i));
            end
            
            % Compute dX2/dx = 2*e'*W*dy/dx = 2*e'*W*C
            val = vec(2*e'*(V(ind,ind)\C1));
        else
            val = zeros(nx, 1);
        end
    end

    function val = d2Gdx2(t,sol)
        nx = size(sol.C1,2);

        % Find all data points that have a time that matches t
        ind = find(t == timelist);
        tlength = length(ind);
        if tlength > 0
            % compute d2X2/dx2 = 2*C1'*W*C1
            C1 = zeros(tlength,nx);
            for i = 1:tlength
                C1(i,:) = sol.C1(outputlist(ind(i)),:);
            end
            val = 2*C1*(V(ind,ind)\C1);
        else
            val = zeros(nx, nx);
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Information theory %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probability density function
    function p = p(sol)
        % p = tau^(-n/2) * prod(sigma)^-1 * exp(-1/2 * sum(((ybar-yhat)/sigma)^2)
        chi2 = G(sol);
        p = (2*pi)^(-n/2) * prod(sigma).^-1 * exp(-1/2 * chi2);
    end

%% Log likelihood
    function val = logp(sol)
        % logp = -n/2 * log(tau) + -sum(log(sigma)) + -1/2 * sum(((ybar-yhat)/sigma)^2)
        chi2 = G(sol);
        val = -n/2 * log(2*pi) + -sum(log(sigma)) + -1/2 * chi2;
    end

%% Fisher information matrix
    function val = F(sol)
        T = [sol.k(sol.UseParams); sol.s(sol.UseSeeds); sol.q(sol.UseControls)];
        nT = numel(T);
        nx = size(sol.C1,2);
        
        % Evaluate the ODE solver structure
        dxdTStart = nx+1;
        dxdTEnd   = nx+nx*nT;
        if isempty(sol.idata)
            % The solution is already discretized
            if numel(sol.x) == numel(discreteTimes)
                dxdT = sol.y(dxdTStart:dxdTEnd,:); % xT_t
            else
                ind = lookup(discreteTimes, sol.x);
                dxdT = sol.y(dxdTStart:dxdTEnd,ind); % xT_t
            end
        else
            % The complete solution is provided
            joint = deval(sol, discreteTimes, 1:dxdTEnd); % x+xT_t
            dxdT = joint(dxdTStart:dxdTEnd,:); % xT_t
        end
        
        % Get dydT for every point
        dydT = zeros(n, nT);
        for i = 1:n
            ind = find(discreteTimes == timelist(i), 1);
            idxdT = reshape(dxdT(:,ind), nx,nT); %x_T
            iC = sol.C1(outputlist(i),:);
            
            % dy(i)/dT
            dydT(i,:) = iC * idxdT;
        end
        
        % Compute expected normalized hessian
        val = dydT.' * (V \ dydT);
        val = symmat(val);
    end

%% Normalized Fisher information matrix
    function val = Fn(sol, T)
        T = [sol.k(sol.UseParams); sol.s(sol.UseSeeds); sol.q(sol.UseControls)];
        nT = numel(T);
        nx = size(sol.C1,2);
        
        % Evaluate the ODE solver structure
        dxdTStart = nx+1;
        dxdTEnd   = nx+nx*nT;
        if isempty(sol.idata)
            % The solution is already discretized
            if numel(sol.x) == numel(discreteTimes)
                dxdT = sol.y(dxdTStart:dxdTEnd,:); % xT_t
            else
                ind = lookup(discreteTimes, sol.x);
                dxdT = sol.y(dxdTStart:dxdTEnd,ind); % xT_t
            end
        else
            % The complete solution is provided
            joint = deval(sol, discreteTimes, 1:dxdTEnd); % x+xT_t
            dxdT = joint(dxdTStart:dxdTEnd,:); % xT_t
        end
        
        % Get dydT for every point
        dydT = zeros(n, nT);
        for i = 1:n
            ind = find(discreteTimes == timelist(i), 1);
            idxdT = reshape(dxdT(:,ind), nx,nT); %x_T
            iC = sol.C1(outputlist(i),:);
            
            % dy(i)/dT
            dydT(i,:) = iC * idxdT;
        end
        
        % Construct normalization matrix
        val = spdiags(T,0,nT,nT); % T along the diagonal

        % Compute expected normalized hessian
        val = dydT*val;
        val = val.' * (V \ val);
        val = symmat(val);
    end

%% AddData
    function objNew = AddData(sol)
        nx = size(sol.C1,1);
        
        % Evaluate the ODE solver structure
        if isempty(sol.idata)
            % The solution is already discretized
            x = sol.y(1:nx, lookup(discreteTimes, sol.x));
            u = sol.u(:,lookup(discreteTimes, sol.x));
        else
            % The complete solution is provided
            x = deval(sol, discreteTimes, 1:nx); % x_t
            u = sol.u(discreteTimes);
        end
        
        % New measurements
        newMeasurements = zeros(n,1);
        for i = 1:n
            t = discreteTimes == timelist(i);
            newMeasurements(i) = sol.C1(outputlist(i),:) * x(:,t) + sol.C2(outputlist(i),:) * u(:,t) + sol.c(outputlist(i),:);
            newMeasurements(i) = newMeasurements(i) + randn*sd(timelist(i),outputlist(i),newMeasurements(i));
        end
        
        newMeasurements(newMeasurements < 0) = 0;
        
        objNew = objectiveWeightedSumOfSquaresFixedToData(outputlist, timelist, sd, newMeasurements, name);
    end
end
