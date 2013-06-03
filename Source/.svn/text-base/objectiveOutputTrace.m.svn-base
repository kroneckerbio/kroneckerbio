function obj = objectiveOutputTrace(m, yInd, yTarget, ySd, name)
% obj = objectiveTargetFunction(m, yInd, yTarget, ySd, name)

%% Work-up
% Constants
ny = m.ny;

% Gut m
temp = m;
clear m
m.ny = temp.ny;
m.C1 = temp.C1;
m.C2 = temp.C2;
m.c  = temp.c;
clear temp

%% Process yInd
if isempty(yInd)
    yInd = true(ny,1);
elseif isnumeric(yInd)
    assert(all(floor(yInd) == yInd) && all(yInd >= 1), 'KroneckerBio:objectiveOutputTarget:InvalidLinearIndex', 'yInd is an invalid linear index')
    assert(all(yInd <= ny), 'KroneckerBio:objectiveOutputTarget:LinearIndexOutOfRange', 'If yInd is provided as a linear index, no element can not be larger than m.ny')
    temp = false(ny,1);
    temp(yInd) = true;
    yInd = temp;
elseif islogical(yInd)
    assert(numel(yInd) <= ny, 'KroneckerBio:objectiveOutputTarget:InvalidLogicalLength', 'If yInd is provided as a logical index, numel(yInd) cannot be larger than m.ny')
    yInd = vec(yInd);
    yInd = [yInd; false(ny-numel(yInd),1)];
else
    error('KroneckerBio:objectiveTargetFunction:IndexType', 'yInd must be a linear or logical index into vectors of length m.ny')
end

nTar = sum(yInd);

%% Process yTarget
% Numeric
if isnumeric(yTarget)
    value = vec(yTarget);
    assert(numel(yTarget) == nTar, 'KroneckerBio:objectiveOutputTarget:TargetWrongLength', 'Numeric targets for yTarget must match the number of indexes yInd')
    yTarget = @(t)repmat(value, 1,numel(t));
end

assert(isa(yTarget, 'function_handle') && nargin(yTarget) == 1, 'KroneckerBio:objectiveOutputTarget:TargetType', 'Target function yTarget must be a vector or a function handle @(t) returning a vector')

%% Process ySd
if isempty(ySd)
    ySd = @(t)ones(nTar,numel(t));
elseif isnumeric(ySd)
    assert(numel(ySd) == nTar || numel(ySd) == 1, 'KroneckerBio:objectiveOutputTarget:WeightWrongLength', 'The numel(ySds) must equal 1 or match the number of indexes yInd')
    value = vec(ySd);
    if numel(value) == 1
        value = repmat(value, nTar,1);
    end
    ySd = @(t)repmat(value, 1,numel(t));
end

assert(isa(ySd, 'function_handle') && nargin(ySd) == 1, 'KroneckerBio:objectiveOutputTarget:WeightType', 'Weight function ySd must be a vector, a scalar, or a function handle @(t) returning a vector')

%% Process name
if isempty(name)
    name = 'Unnamed Objective';
end

%% Contruct objective function
obj.Type = 'Objective.Data.ChiSquare';
obj.Name = name;
obj.DiscreteTimes = zeros(0,1);
obj.Continuous = true;
obj.Complex = false;
obj.Linked = 0;

obj.g = @g;
obj.dgdx = @dgdx;
obj.d2gdx2 = @d2gdx2;

obj.Update = @Update;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function val = g(t, x ,u)
        % g = (y - ytar) *{y.y} (sigma^-1 *{sxy} (y - ytar))
        % Value of y(t)
        val = m.C1(yInd,:) * x + m.C2(yInd,:) * u + m.c(yInd,:);
        
        % Weighted difference
        val = (val - yTarget(t)) ./ ySd(t);
        
        % Sum of squares
        val = val.' * val;
    end

    function val = dgdx(t, x ,u)
        % dgdx = 2 * dydx * (sigma^-1 *{sxy} (y - ytar))
        % Value of y(t)
        val = m.C1(yInd,:) * x + m.C2(yInd,:) * u + m.c(yInd,:);
        
        % Weighted difference
        val = (val - yTarget(t)) ./ ySd(t);
        
        % dydx = C1
        val = 2 * (m.C1(yInd,:).' * val);
        
        val = vec(val);
    end

    function val = d2gdx2(t, x ,u)
        % d2gdx2 = 2 * dydx *{y.y} (sigma^-1 *{sxy} dydx)
        val = 2 * m.C1(yInd,:).' * bsxfun(@rdivide, m.C1(yInd,:), ySd(t));
    end

%% Update
    function objNew = Update(m, con, UseParams, UseICs, UseControls)
        objNew = objectiveOutputTrace(m, yInd, yTarget, ySd, name);
    end
end