function obs = observationLogWeightedSumOfSquares(outputlist, timelist, sd, name)
%observationLogWeightedSumOfSquares Create an observation scheme for a
%   list of points, whose uncertainty is a log normal distribution with a
%   central tendency of the true value and a standard deviation that is an
%   arbitrary function of the true value.
%
%   obs = observationLogWeightedSumOfSquares(outputlist, timelist, sd, 
%                                            name)
%
%   Inputs
%   outputlist: [ positive integer vector n ]
%       The indexes of the outputs for each measurement.
%   timelist: [ nonnegative vector n ]
%       The times for each measurement
%   sd: [ handle @(t,i,y) returns positive ]
%       The standard deviation function that takes the time, the index, and
%       the true value and returns the standard deviation of the
%       measurement. Because there is ambiguity in the standard deviation
%       of a log normal distribution, this function specifically expects a
%       return value equal to (sigma_logy * y), where sigma_logy is the
%       standard deviation of log(y) and y is normal output value.
%   name: [ string ]
%       An arbitrary name for the observation scheme
%
%   Outputs
%   obs: [ observation scheme structure ]

% (c) 2017 David R Hagen
% This work is released under the MIT license.

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

obs.Type = 'Observation.Data.LogWeightedSumOfSquares';
obs.Name = name;
obs.Complex = false;

obs.tF = max([0;timelist]);
obs.DiscreteTimes = discrete_times;

obs.Simulation = @simulation;
obs.Objective  = @objective;

obs.F = @(sol)F(outputlist, timelist, discrete_times, sd, sol);

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        y_all = int.y;
        
        ybar = zeros(n,1);
        yhat = zeros(n,1);
        for i = 1:n
            ind = discrete_times == timelist(i);
            ybar(i) = y_all(outputlist(i),ind);
            yhat(i) = exp(log(ybar(i)) + randn * sd(timelist(i), outputlist(i), ybar(i)) / ybar(i)); % linear -> log -> linear
        end
        
        sim.Type = 'Simulation.Data.LogWeightedSumOfSquares';
        sim.Name = name;
        sim.int = int;
        sim.outputlist = outputlist;
        sim.timelist = timelist;
        sim.true_measurements = ybar;
        sim.measurements = yhat;
        sim.sd = sd;
    end

    function obj = objective(measurements)
        measurements = vec(measurements);
        assert(numel(measurements) == n , 'KroneckerBio:observationWeightedSumOfSquares:measurements', 'Input "measurements" must be a vector length of "outputlist"')
        obj = objectiveLogWeightedSumOfSquares(outputlist, timelist, measurements, sd, name);
    end
end

function obj = objectiveLogWeightedSumOfSquares(outputlist, timelist, measurements, sd, name)
% Log measurements
log_measurements = log(measurements);

% Find unique timelist
n = numel(outputlist);
discrete_times = row(unique(timelist));

% Inherit observation
obj = observationLinearWeightedSumOfSquares(outputlist, timelist, sd, name);

obj.Type = 'Objective.Data.LogWeightedSumOfSquares';
obj.Continuous = false;

obj.G = @G;
obj.dGdy = @dGdy;
obj.d2Gdy2 = @d2Gdy2;

obj.p = @p;
obj.logp = @logp;
obj.F = @(sol)F(outputlist, timelist, discrete_times, sd, sol);

obj = pastestruct(objectiveZero(), obj);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val, discrete] = G(int)
        % Evaluate solution
        [logybar, logsigma] = evaluate_sol(outputlist, timelist, discrete_times, sd, int);
        loge = logybar - log_measurements;
        
        % Goal function
        val = sum(2 * log(logsigma) + (loge./logsigma).^2);
        
        % Return discrete times as well
        discrete = discrete_times;
    end

    function val = dGdy(t, int)
        % dGdy = 2 * sigma^-1 * dsigmady + 2 * (y-yhat) * (sigma - (y-yhat)*dsigmady) * sigma^-3
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
            log_measurements_t = log_measurements(ind_t);
            
            ybar_t = zeros(n_current,1);
            logybar_t = zeros(n_current,1);
            logsigma_t = zeros(n_current,1);
            dlogsigmadlogy_t = zeros(n_current,1);
            for i = 1:n_current
                ybar_t(i) = yt(outputlist_t(i));
                logybar_t(i) = log(ybar_t(i));
                [sigma_t_i, dsigmady_t_i] = sd(timelist_t(i), outputlist_t(i), ybar_t(i));
                logsigma_t(i) = sigma_t_i ./ ybar_t(i);
                dlogsigmadlogy_t(i) = dsigmady_t_i - logsigma_t(i);
            end
            
            % Gradient value
            e = logybar_t - log_measurements_t; % Y_
            dGdlogybar_t = 2 ./ logsigma_t .* dlogsigmadlogy_t + 2 .* e .* (logsigma_t - e.*dlogsigmadlogy_t) ./ logsigma_t.^3; % Y_
            val = accumarray(outputlist_t, dGdlogybar_t ./ ybar_t, [ny,1]); % y_ % sum the entries associated with the same output
        else
            val = zeros(ny, 1);
        end
    end

    % TODO: there is something wrong with this  derivative.
    function val = d2Gdy2(t, int)
        error('observationLogWeightedSumOfSquares.d2Gdy2 is inaccurate and needs fixed')
        % d2Gdy2dy1 = 2 * sigma^-2 - 4*e*sigma^-3*dsigmady1 - 4*e*sigma^-3*dsigmady2 
        %             + (6*e^2*sigma^-4 - 2*sigma^-2)*dsigmady1*dsigmady2 
        %             + (2*sigma^-1 - 2*e^2*sigma^-3)*d2sigmady1dy2
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
            log_measurements_t = log_measurements(ind_t);
            
            ybar_t = zeros(n_current,1);
            logybar_t = zeros(n_current,1);
            logsigma_t = zeros(n_current,1);
            dlogsigmadlogy_t = zeros(n_current,1);
            d2logsigmadlogy2_t = zeros(n_current,1);
            for i = 1:n_current
                ybar_t(i) = yt(outputlist_t(i));
                logybar_t(i) = log(ybar_t(i));
                [sigma_t_i, dsigmady_t_i, d2sigmady2_t_i] = sd(timelist_t(i), outputlist_t(i), ybar_t(i));
                logsigma_t(i) = sigma_t_i ./ ybar_t(i);
                dlogsigmadlogy_t(i) = dsigmady_t_i - logsigma_t(i);
                d2logsigmadlogy2_t(i) = d2sigmady2_t_i .* ybar_t(i) - 2 * dlogsigmadlogy_t(i); % Not sure if this factor of 2 should be here
            end
            
            % Curvature value
            loge = logybar_t - log_measurements_t;
            d2Gdlogybar2_t = 2 ./ logsigma_t.^2 - 8 .* loge ./ logsigma_t.^3 .* dlogsigmadlogy_t + ...
                (6 .* loge.^2 ./ logsigma_t.^4 - 2 ./ logsigma_t.^2) .* dlogsigmadlogy_t.^2 + ...
                (2 ./ logsigma_t - 2 .* loge.^2 ./ logsigma_t.^3) .* d2logsigmadlogy2_t; % Y_
            val = accumarray(outputlist_t, d2Gdlogybar2_t ./ ybar_t ./ ybar_t, [ny,1]); % y_ % sum the entries associated with the same output
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
        [logybar, logsigma] = evaluate_sol(outputlist, timelist, discrete_times, sd, sol);
        loge = logybar - log_measurements;
        
        val = (2*pi)^(-n/2) * prod(logsigma).^-1 * exp(-1/2 * sum((loge ./ logsigma).^2));
    end

%% Log likelihood
    function val = logp(sol)
        % logp = -n/2 * log(tau) + -sum(log(sigma)) + -1/2 * sum(((ybar-yhat)/sigma)^2)

        % Evaluate solution
        [logybar, logsigma] = evaluate_sol(outputlist, timelist, discrete_times, sd, sol);
        loge = logybar - log_measurements;
        
        val = -n/2 * log(2*pi) + -sum(log(logsigma)) + -1/2 * sum((loge ./ logsigma).^2);
    end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Helper functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [logybar, logsigma] = evaluate_sol(outputlist, timelist, discrete_times, sd, int)
n = numel(outputlist);
y_all = int.y;

ybar = zeros(n,1);
sigma = zeros(n,1);
for i = 1:n
    ind = discrete_times == timelist(i);
    ybar(i) = y_all(outputlist(i),ind);
    sigma(i) = sd(timelist(i), outputlist(i), ybar(i));
end

% Convert to log space
logybar = log(ybar);
logsigma = sigma ./ ybar;
end

function [dlogydT, logsigma, dlogsigmadlogy] = evaluate_grad(outputlist, timelist, discrete_times, sd, int)
n = numel(outputlist);
nT = int.nT;
ny = int.ny;

y_all = int.y;
dydT_all = int.dydT;

% Get dydT for every point
ybar = zeros(n,1);
dydT = zeros(n, nT);
sigma = zeros(n,1);
dsigmady = zeros(n,1);
for i = 1:n
    ind = find(discrete_times == timelist(i), 1);
    t_i = timelist(i);
    
    % Compute this y
    ybar(i) = y_all(outputlist(i),ind);
    
    % Compute dy/dT
    dydT_temp = reshape(dydT_all(:,ind), ny, nT);
    dydT(i,:) = dydT_temp(outputlist(i),:);
    
    % Compute expected V matrix
    [sigma(i), dsigmady(i)] = sd(t_i, outputlist(i), ybar(i));
end

dlogydT = bsxfun(@rdivide, dydT, y_all);
logsigma = sigma ./ ybar;
dlogsigmadlogy = dsigmady - logsigma;
end

%% Fisher information matrix
function val = F(outputlist, timelist, discrete_times, sd, sol)
n = numel(outputlist);

% Evaluate solution
[dlogydT, logsigma, dsigmadlogy] = evaluate_grad(outputlist, timelist, discrete_times, sd, sol);
dlogsigmadT = diag(dsigmadlogy)*dlogydT;
dCovdlogsigma = 2*diag(logsigma);
dCovdT = dCovdlogsigma*dlogsigmadT;
nT = size(dlogydT,2);

% Assemble variance matrix
V = spdiags(logsigma.^2,0,n,n);

% Fisher information matrix
covarTerm = zeros(nT);
diagCovInv = diag(logsigma.^-2);
for m = 1:nT
    dCovdT_m = diag(dCovdT(:,m));
    for n = 1:nT
        dCovdT_n = diag(dCovdT(:,n));
        covarTerm(m,n) = 1/2*trace(diagCovInv*dCovdT_m*diagCovInv*dCovdT_n);
    end
end
val = dlogydT.' * (V \ dlogydT) + covarTerm;
val = symmat(val);
end
