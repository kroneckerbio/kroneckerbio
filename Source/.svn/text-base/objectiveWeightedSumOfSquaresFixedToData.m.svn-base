function obj = objectiveWeightedSumOfSquaresFixedToData(m, outputlist, timelist, sd, nonNegMeasurements, measurements, perfect, name)
%obj = objectiveWeightedSumOfSquaresFixedToData(m, outputlist, timelist, sd,
%nonNegMeasurements, measurements, perfect)
%returns X2, chi-square, the generalized least squares statistic

% sd: @(t,yInd,yVal) returns standard error given t, index of output, and
% value of output

% Clean up inputs
if nargin < 8
    name = [];
    if nargin < 7
        perfect = [];
        if nargin < 6
            measurements = [];
            if nargin < 5
                nonNegMeasurements = [];
                if nargin < 4
                    sd = [];
                    if nargin < 3
                        timelist = [];
                        if nargin < 2
                            outputlist = [];
                            if nargin < 1
                                m = [];
                            end
                        end
                    end
                end
            end
        end
    end
end

% Special case: return empty structure array if inputs are numeric
if isnumeric(m)
    obj = emptystruct(m, 'Type', 'Name', 'DiscreteTimes', 'Continuous', 'Complex', 'G', 'dGdx', 'd2Gdx2', 'F', 'Fn', 'dFndT', 'p', 'pvalue', 'n', 'AddData', 'AddExpectedData', 'QuantifyPrediction', 'Update');
    return
end

% Check inputs
assert(isvector(outputlist), 'KroneckerBio:constructObjectiveChiSquare:outputlist', 'Input "outputlist" must be a vector.')
n = length(outputlist);
assert(isvector(timelist) && length(timelist) == n, 'KroneckerBio:constructObjectiveChiSquare:timelist', 'Input "timelist" must be a vector length of "outputlist".')
assert(numel(measurements) == n || isempty(measurements), 'KroneckerBio:constructObjectiveChiSquare:measurements', 'Input "measurements" must be a vector length of "outputlist" or empty.')

% Gut m
temp = m;
clear m
m.nx = temp.nx;
m.nk = temp.nk;
m.c = temp.c;
clear temp

% Constants
nx = m.nx;
nk = m.nk;

% Find unique times
discreteTimes = vec(unique(timelist)).';

% Compute covariance matrix
if ~isempty(measurements)
    V = zeros(n,1);
    for j = 1:n
        V(j) = sd(timelist(j), outputlist(j), measurements(j))^2;
    end
    V = sparse(1:n, 1:n, V, n,n);
end

% Default is nonNegMeasurements
if isempty(nonNegMeasurements)
    nonNegMeasurements = true;
end

% Default is real data
if isempty(perfect)
    perfect = false;
end

if isempty(name)
    name = 'ChiSquare';
end

obj.Type = 'Objective.Data.ChiSquare';
obj.Name = name;
obj.DiscreteTimes = discreteTimes;
obj.Continuous = false;
obj.Complex = false;
obj.Linked = 0;

obj.G = @G;
obj.dGdx = @dGdx;
obj.d2Gdx2 = @d2Gdx2;

obj.F = @F;
obj.Fn = @Fn;
obj.dFndT = @dFndT;
obj.p = @p;
obj.pvalue = @pvalue;

obj.n = n;
obj.AddData = @AddData;
obj.AddExpectedData = @AddExpectedData;
obj.QuantifyPrediction = @QuantifyPrediction;

obj.Update = @update;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [val discrete] = G(sol)
        if ~perfect
            % Evaluate the ODE solver structure
            if isempty(sol.idata)
                % The solution is already discretized
                if numel(sol.x) == numel(discreteTimes)
                    x = sol.y(1:nx,:); % x_t
                    u = sol.u;
                else
                    ind = lookup(discreteTimes, sol.x);
                    x = sol.y(1:nx,ind);
                    u = sol.u(:,ind);
                end
            else
                % The complete solution is provided
                x = deval(sol, discreteTimes, 1:nx); % x_t
                u = sol.u(discreteTimes);
            end
            
            % Get e for every point evaluated
            e = zeros(n,1);
            for i = 1:n
                t = discreteTimes == timelist(i);
                e(i,1) = sol.C1(outputlist(i),:) * x(:,t) + sol.C2(outputlist(i),:) * u(:,t) + sol.c(outputlist(i),:) - measurements(i);
            end
            
            % Sort by absolute value
            [unused I] = sort(abs(e.'*sqrt(V)));
            e = e(I);
            Vsort = diag(V); % Extract diagonal
            Vsort = spdiags(Vsort(I),0,n,n); % Replace diagonal with sorted version
            
            % Compute chi-square
            val = e'*(Vsort\e);
            
            if nargout >= 2
                discrete = discreteTimes;
            end
        else
            % Perfect data requires an expected goal value
            val = EG(sol);
            discrete = discreteTimes;
        end
    end

    function val = dGdx(t,sol)
        %Find all data points that have a time that matches t
        ind = find(t == timelist);
        tlength = length(ind);
        if tlength > 0
            % Evaluate the ODE solver structure
            if isempty(sol.idata)
                % The solution is already discretized
                % Get the solution at the matching timepoint
                tInd = sol.x == t;
                x = sol.y(1:nx, tInd);
                u = sol.u(:,tInd);
            else
                % The complete solution is provided
                x = deval(sol, t, 1:nx);
                u = sol.u(t);
            end
            
            % Compute e for each datapoint that has a matching time to t
            e = zeros(tlength,1);
            C1 = zeros(tlength,nx);
            for i = 1:tlength
                C1(i,:) = sol.C1(outputlist(ind(i)),:);
                e(i,1) = C1(i,:)*x + sol.C2(outputlist(ind(i)),:) * u + sol.c(outputlist(ind(i)),:) - measurements(ind(i));
            end
            
            % Sort by absolute value
            Vsort = V(ind,ind);
            [unused I] = sort(abs(e.'*sqrt(Vsort)));
            e = e(I);
            C1 = C1(I,:);
            Vsort = diag(Vsort); % Extract diagonal
            Vsort = spdiags(Vsort(I),0,tlength,tlength); % Replace diagonal with sorted version
            
            % Compute dX2/dx = 2*e'*W*dy/dx = 2*e'*W*C
            val = vec(2*e'*(Vsort\C1));
        else
            val = zeros(nx, 1);
        end
    end

    function val = d2Gdx2(t,sol)
        %Find all data points that have a time that matches t
        ind = find(t == timelist);
        tlength = length(ind);
        if tlength > 0
            %compute d2X2/dx2 = 2*C1'*W*C1
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
%% Fisher information matrix
    function val = F(sol)
        nT = (size(sol.y, 1) - nx) / nx;
        
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
        nT = (size(sol.y, 1) - nx) / nx;
        
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

    function val = dFndT(sol, T)
        nT = numel(T);
        
        % Evaluate the ODE solver structure
        dxdTStart   = nx+1;
        dxdTEnd     = nx+nx*nT;
        d2xdT2Start = nx+nx*nT+1;
        d2xdT2End   = nx+nx*nT+nx*nT*nT;
        if isempty(sol.idata)
            % The solution is already discretized
            if numel(sol.x) == numel(discreteTimes)
                dxdT = sol.y(dxdTStart:dxdTEnd,:); % xT_t
                d2xdT2 = sol.y(d2xdT2Start:d2xdT2End,:); % xTT_t
            else
                ind = lookup(discreteTimes, sol.x);
                dxdT = sol.y(dxdTStart:dxdTEnd,ind); % xT_t
                d2xdT2 = sol.y(d2xdT2Start:d2xdT2End,ind); % xTT_t
            end
        else
            % The complete solution is provided
            joint = deval(sol, discreteTimes, 1:d2xdT2End); % x+xT_t
            dxdT = joint(dxdTStart:dxdTEnd,:); % xT_t
            d2xdT2 = joint(d2xdT2Start:d2xdT2End,:); % xTT_t
        end
        
        % Get dydT for every point
        dydT = zeros(n, nT);
        d2ydT2 = zeros(n, nT*nT);
        for i = 1:n
            ind = find(discreteTimes == timelist(i), 1);
            idxdT = reshape(dxdT(:,ind), nx,nT); %x_T
            id2xdT2 = reshape(d2xdT2(:,ind), nx,nT*nT); % x_TT
            iC = sol.C1(outputlist(i),:);
            
            % dy(i)/dT
            dydT(i,:) = iC * idxdT;
            
            % d2y(i)/dT2
            d2ydT2(i,:) = iC * id2xdT2;
        end
        
        dydTlog = bsxfun(@times, dydT, T.'); % T is normalized
        
        d2ydT2 = reshape(d2ydT2, n*nT,nT);
        d2ydT2 = spermute132(bsxfun(@times, d2ydT2, T.'), [n,nT,nT], [n,nT*nT]); % Inner T is normalized
        
        dydTdiag = sparse(repmat(vec(1:n), 1,nT), repmat((0:nT+1:nT*nT-1)+1,n,1), dydT);
        
        val = dydTlog.' * (V \ d2ydT2) + dydTlog.' * (V \ reshape(dydTdiag, n,nT*nT));
        val = val + spermute213(val, [nT,nT,nT], [nT,nT*nT]);
        val = reshape(val, nT*nT,nT);
    end

%% Probability density function
    function p = p(sol)
        % This is a centralized probability density function. That is, when
        % chi2 = n, then p = 1. The probability has been scaled to
        % accomodate this and is necessary to store the probability as a
        % floating point number. The relative probabilities between
        % solutions using the same objective function (or even objective
        % functions with the same number of data points) remains correct.
        % p is a direct consequence of chi-square
        
        chi2 = G(sol);
        %p = (2*pi)^(-n/2) * det(V)^(-1/2) * exp(-1/2*chi2); % non-centralized
        p = exp(-1/2 * (chi2-n)); % centralized
    end

%% P-Value
    function pval = pvalue(sol, obj)
        nCon = numel(sol);
        nObj = size(obj,1);
        chi2 = 0;
        ntot = 0;
        
        % Sum chi-square and n across all data
        for iCon = 1:nCon
            for iObj = 1:nObj
                chi2 = chi2 + obj(iObj,iCon).G(sol(iCon));
                ntot = ntot + obj(iObj,iCon).n;
            end
        end
        
        pval = chi2pvalue(chi2, ntot);
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
        
        if nonNegMeasurements
            newMeasurements(newMeasurements < 0) = 0;
        end
        
        objNew = constructObjectiveChiSquare(m, outputlist, timelist, sd, nonNegMeasurements, newMeasurements, false);
    end

%% AddExpectedData
    function objNew = AddExpectedData(sol)
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
            tInd = discreteTimes == timelist(i);
            newMeasurements(i) = sol.C1(outputlist(i),:) * x(:,tInd) + sol.C2(outputlist(i),:) * u(:,tInd) + sol.c(outputlist(i),:);
        end
        
        % Build new objective function with perfect data
        objNew = constructObjectiveChiSquare(m, outputlist, timelist, sd, nonNegMeasurements, newMeasurements, true);
    end

%% QuantifyPrediction
    function result = QuantifyPrediction(sol1, sol2, type)
        % Special case: return empty structure array if inputs are numeric
        if isnumeric(sol1)
            result = zeros(sol1);
            return
        end
        
        % Check type
        if isempty(type)
            type = 'diffscalemaxmean';
        end
        
        % Evaluate the ODE solver structure
        if isempty(sol1.idata)
            % The solution is already discretized
            if all(sol1.x == discreteTimes)
                x1 = sol1.y(1:nx,:); % x_t
                x2 = sol2.y(1:nx,:); % x_t
            else
                error('Need to implement multiple obj functions here')
            end
        else
            % The complete solution is provided
            x1 = deval(sol1, discreteTimes, 1:nx); % x_t
            x2 = deval(sol2, discreteTimes, 1:nx); % x_t
        end
        
        % Convert to outputs
        y1 = zeros(n,1);
        y2 = zeros(n,1);
        for i = 1:n
            y1(i,1) = m.c(outputlist(i),:)*x1(:,discreteTimes == timelist(i));
            y2(i,1) = m.c(outputlist(i),:)*x2(:,discreteTimes == timelist(i));
        end
        
        switch lower(type)
            case 'diffscalemaxmean'
                max1 = max(x1,[],2).^-1;
                max1 = spdiags(max1, 0, nx,nx);
                
                % Scale species by the maximum
                x1 = max1 * x1;
                x2 = max1 * x2;
                
                % Convert to outputs
                y1 = zeros(n,1);
                y2 = zeros(n,1);
                for i = 1:n
                    y1(i,1) = m.c(outputlist(i),:)*x1(:,discreteTimes == timelist(i));
                    y2(i,1) = m.c(outputlist(i),:)*x2(:,discreteTimes == timelist(i));
                end
                
                % Difference
                result = mean(abs(y1 - y2), 1);
                
            case 'diffscalemaxmax'
                max1 = max(x1,[],2).^-1;
                max1 = spdiags(max1, 0, nx,nx);
                
                % Scale species by the maximum
                x1 = max1 * x1;
                x2 = max1 * x2;
                
                % NaN should be zero
                x1(isnan(x1)) = 0;
                x2(isnan(x2)) = 0;
                
                % Convert to outputs
                y1 = zeros(n,1);
                y2 = zeros(n,1);
                for i = 1:n
                    y1(i,1) = m.c(outputlist(i),:)*x1(:,discreteTimes == timelist(i));
                    y2(i,1) = m.c(outputlist(i),:)*x2(:,discreteTimes == timelist(i));
                end
                
                % Difference
                result = max(abs(y1 - y2), [], 1);
                
            case 'diffscalemaxtruemax'
                max1 = max(x1,[],2).^-1;
                max1 = spdiags(max1, 0, nx,nx);
                
                % Scale species by the maximum
                x1 = max1 * x1;
                x2 = max1 * x2;
                
                % NaN should be zero
                x1(isnan(x1)) = 0;
                x2(isnan(x2)) = 0;
                
                % Convert to outputs
                y1 = zeros(n,1);
                y2 = zeros(n,1);
                for i = 1:n
                    y1(i,1) = m.c(outputlist(i),:)*x1(:,discreteTimes == timelist(i));
                    y2(i,1) = m.c(outputlist(i),:)*x2(:,discreteTimes == timelist(i));
                end
                
                % Difference
                result = max(abs(y1 - y2), [], 1);
                
        end % switch
    end

%% Update
    function objNew = update(m, con, UseParams, UseICs, UseControls)
        objNew = objectiveWeightedSumOfSquaresFixedToData(m, outputlist, timelist, sd, nonNegMeasurements, measurements, perfect, name);
    end

%% Expected goal function
% The expected value of normalized chi-square given the current
% measurements as the bias. If measurements is unset, there is assumed to
% be no bias.
    function val = EG(sol)
        % Evaluate the ODE solver structure
        if isempty(sol.idata)
            % The solution is already discretized
            if numel(sol.x) == numel(discreteTimes)
                x = sol.y(1:nx,:); % x_t
                u = sol.u;
            else
                ind = lookup(discreteTimes, sol.x);
                x = sol.y(1:nx,ind);
                u = sol.u(:,ind);
            end
        else
            % The complete solution is provided
            x = deval(sol, discreteTimes, 1:nx); % x_t
            u = sol.u(discreteTimes);
        end
        
        % Get bias d for every point evaluated
        d = zeros(n,1);
        for i = 1:n
            t = discreteTimes == timelist(i);
            d(i,1) = sol.C1(outputlist(i),:) * x(:,t) + sol.C2(outputlist(i),:) * u(:,t) + sol.c(outputlist(i),:) - measurements(i);
        end
        
        % Sort by absolute value
        [unused I] = sort(abs(d.'*sqrt(V)));
        d = d(I);
        Vsort = diag(V); % Extract diagonal
        Vsort = spdiags(Vsort(I),0,n,n); % Replace diagonal with sorted version

        % Compute expected chi-square
        val = n + d.' * (Vsort \ d);
    end
end