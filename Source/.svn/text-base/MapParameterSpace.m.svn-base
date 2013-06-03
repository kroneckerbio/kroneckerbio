function varargout = MapParameterSpace(m, con, obj, opts)
% TODO: allow heat map for march method
% TODO: allow marching to march in parameter directions and enumeration to
% evaluate along eigenvectors
% TODO: fix subplot --> axes for march method
% TODO: check that zpvalue is the appropriate pvalue function for
% chi-square (or anything really)

%% Work-up
% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:MapParameterSpace:TooFewInputs', 'MapParameterSpace requires at least 3 input arguments')
if nargin < 4
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:MapParameterSpace:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseICs         = [];
defaultOpts.UseControls    = [];

defaultOpts.ObjWeights     = ones(size(obj));

defaultOpts.Normalized     = true;
defaultOpts.UseAdjoint     = true;

defaultOpts.SampleCount     = 3;
defaultOpts.Pairs           = 'all';
defaultOpts.ToPlot          = false;
defaultOpts.PlotProgress    = false;
defaultOpts.PlotType        = 'contour';
% Marching
defaultOpts.PValueFun       = @zpvalue;
defaultOpts.PValueTargets   = 0.32;
defaultOpts.Normalized      = true;
defaultOpts.MaxLength       = 1000;
% Enumeration
defaultOpts.EnumSource      = 'real';
defaultOpts.LowerBound      = m.k * 0.1;
defaultOpts.UpperBound      = m.k * 10;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nk = m.nk;
nCon = numel(con);
nObj = size(obj,1);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseICs is a logical matrix
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseControls nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon, m.nq, cat(1,con.nq));

nT = nTk + nTx + nTq;

opts.LowerBound = fixBounds(opts.LowerBound, opts.UseParams, opts.UseICs, opts.UseModelICs);
opts.UpperBound = fixBounds(opts.UpperBound, opts.UseParams, opts.UseICs, opts.UseModelICs);

nSam = opts.SampleCount;
nCut = length(opts.PValueTargets);

%% Create all sets of pairs
if ischar(opts.Pairs)
    if strcmpi(opts.Pairs, 'all')
        pairIndexes = combinator(nT, 2, 'c', 'n');
    elseif strcmpi(opts.Pairs, 'adjacent')
        pairIndexes = [vec(1:(nT-1)), vec(2:nT)];
    elseif strcmpi(opts.Pairs, 'first')
        pairIndexes = [1,2];
    else
        error('KroneckerBio:MapParameterSpace:PairingType', 'Pairs must be ''all'', ''adjacent'', ''first'', or numeric')
    end
else%isnumeric
    pairIndexes = opts.Pairs;
end
nPair = size(pairIndexes, 1);

if strcmpi(opts.MapMethod, 'march')
%% Marching method
    %% Find chi-square cutoffs
    pvalEvalFun = opts.PValueFun;
    
    % All p-value cutoffs must be sorted into a column and sorted descending
    pvalCutoffs = opts.PValueTargets ./ 2;
    pvalCutoffs = sort(vec(pvalCutoffs), 1, 'descend');
    
    % Solve for zeros to get chi2 p-value cutoffs
    goalCutoffs = zeros(nCut,1);
    for i = 1:nCut
        goalCutoffs(i) = fzero(@(x)(pvalEvalFun(x) - pvalCutoffs(i)), 0);
    end
    
    %% Compute the fisher information matrix
    F = ObjectiveInformation(m, con, obj, opts);
    
    % Fetch eigenvectors of the fisher information matrix
    [V lambda] = eig(F);
    lambda = diag(lambda);
    lambda(lambda < eps) = eps;
    
    %% Find points at cutoffs
    
    % Enumerate sampling points: 0 to 2pi not repeating at 0
    theta = linspace(0, 2*pi*(nSam-1)/nSam, nSam);
    
    % Starting value of objective function
    p = [m.k(opts.UseParams); m.x0(opts.UseICs)];
    G0 = ObjectiveValue(m, con, obj, opts);
    
    % Create storage for data
    vectors = zeros(2, nPair, nSam);
    radii = zeros(nPair, nSam, nCut);
    flags = zeros(nPair, nSam, nCut);
    
    % Cycle through each pair of eigenvectors
    for iPair = 1:nPair
        % Select one pair
        curPair = V(:, pairIndexes(iPair,:));
        
        % Compute scaling factor for vectors
        a = 1 / sqrt(lambda(pairIndexes(iPair,1))); % ellipse axis 1
        b = 1 / sqrt(lambda(pairIndexes(iPair,2))); % ellipse axis 2
        
        if opts.PlotProgress
            iFig = figure();
        end
        
        % Cycle through each sampling point
        for iSam = 1:nSam
            % ***Create new unit vector***
            % Compute direction on unit circle
            x = cos(theta(iSam));
            y = sin(theta(iSam));
            
            % Create vector pointing that point on a ellipse
            u = [a; b] .* [x; y];
            % Scale back to unit vector
            u = u ./ norm(u);
            
            % Store unit vector for plotting
            vectors(:, iPair, iSam) = u;
            
            % Multiply into eigenvectors
            v = curPair * u;
            
            % ***Find limit to radius***
            % Define boundary for l to prevent p from going to zero
            if opts.Normalized
                limits = -1 ./ v;           % Solve for length needed to make parameter go to zero
            else
                limits = -p ./ v;
            end
            limits(limits < 0) = NaN;   % Destroy the negative ones
            closestLimit = min(limits); % The smallest limit is our boundary for l
            
            % Compute uncertainty in combo direction of eigenvalue pair
            d = sqrt( (a*x)^2 + (b*y)^2 );
            
            % Max distance is a multiple of this uncertainty
            maxLength = d * opts.MaxLength;
            closestLimit = min([closestLimit maxLength]); % Apply to limit
            
            % Find cutoffs
            for iCut = 1:nCut
                % Check if the boundary is sufficient
                Glimit = optFun(closestLimit);
                if Glimit < 0
                    % fzero will never converge
                    l = closestLimit;
                    if closestLimit == maxLength
                        flag = 2; % MaxLength reached
                    else
                        flag = 1; % Zero parameter reached
                    end
                else
                    % fzero is guarenteed to converge
                    l0 = [0 closestLimit];
                    l = fzero(@optFun, l0);
                    flag = 0; % Solution found
                end
                
                % Store length to cutoff
                flags(iPair, iSam, iCut) = flag;
                radii(iPair, iSam, iCut) = l;
                
                % Plot progress
                if opts.PlotProgress
                    figure(iFig)
                    hold on
                    plotx = l * u(1);
                    ploty = l * u(2);
                    plot(plotx, ploty, 'o')
                    drawnow
                end
            end
        end
    end
    
    %% Plot results
    if opts.ToPlot || nargout == 0
        % Plot everything
        for iPair = 1:nPair
            subplot(nT-1, nT-1, (pairIndexes(iPair,1)-1)*(nT-1)+(pairIndexes(iPair,2)-1))
            width = 1 / (nT - 1);
            set(gca, 'ActivePositionTroperty', 'outerposition')
            set(gca, 'OuterPosition', [(pairIndexes(iPair,2)-2)*(width), 1-pairIndexes(iPair,1)*(width), width, width])
            set(gca, 'Position', [(pairIndexes(iPair,2)-2)*(width), 1-pairIndexes(iPair,1)*(width), width, width])
            for iCut = 1:nCut
                % Plot every sample from a particular cut of a pair
                x = radii(iPair,:,iCut) .* reshape(vectors(1,iPair,:), 1,nSam);
                x(end+1) = x(1); % Complete the circle
                y = radii(iPair,:,iCut) .* reshape(vectors(2,iPair,:), 1,nSam);
                y(end+1) = y(1); % Complete the circle
                plot(x,y,'-o')
                hold on
            end
            set(gca, 'XTick', [])
            set(gca, 'YTick', [])
            hold off
        end
    end
    
    %% Work-down
    if nargout >= 1
        % Information matrix
        data.F              = F;
        data.Eigenvalues    = lambda;
        data.Eigenvectors   = V;
        % Independent Variables
        data.Eigenpairs     = pairIndexes;
        data.PValueCutoffs  = pvalCutoffs;
        data.Angles         = theta;
        data.Vectors        = vectors;
        % Length
        data.Radii          = radii;
        data.Flags          = flags;
        
        varargout{1} = data;
    end


else%SearchMethod == 'enum'
%% Enumeration method
    % Use bounds and resolution to compute the enumeration values
    paramEnums = zeros(nT, nSam);
    for iT = 1:nT
        if opts.Normalized
            paramEnums(iT,:) = logspace(log10(opts.LowerBound(iT)), log10(opts.UpperBound(iT)), nSam);
        else
            paramEnums(iT,:) = linspace(opts.LowerBound(iT), opts.UpperBound(iT), nSam);
        end
    end
    
    % Initialize cell array for holding all parameter sets to be tested
    paramSets = cell(nPair, nSam, nSam);
    for iPair = 1:nPair
        % Enumerate the parameter sets that we are going to evaluate
        for iSam = 1:nSam % x-axis
            for jSam = 1:nSam % y-axis
                iEnum = [m.k; m.x0];
                if pairIndexes(iPair,1) <= nTk
                    % Parameter is a rate constant
                    iEnum(useParamsInd(pairIndexes(iPair,1))) = paramEnums(pairIndexes(iPair,1), iSam);
                else
                    % Parameter is an initial condition
                    iEnum(useICsInd(pairIndexes(iPair,1)-nTk)+m.nk) = paramEnums(pairIndexes(iPair,1), iSam);
                end
                if pairIndexes(iPair,2) <= nTk
                    % Parameter is a rate constant
                    iEnum(useParamsInd(pairIndexes(iPair,2))) = paramEnums(pairIndexes(iPair,2), jSam);
                else
                    % Parameter is an initial condition
                    iEnum(useICsInd(pairIndexes(iPair,2)-nTk)+m.nk) = paramEnums(pairIndexes(iPair,2), jSam);
                end
                paramSets{iPair, iSam, jSam} = iEnum;
            end
        end
    end
    
    % Compute the FIM or hessian if needed
    if strcmpi(opts.EnumSource, 'FIM') 
        H = ObjectiveInformation(m, con, obj, opts);
        G0 = ObjectiveValue(m, con, obj, opts);
    elseif strcmpi(opts.EnumSource, 'hessian')
        H = ObjectiveHessian(m, con, obj, opts);
        G0 = ObjectiveValue(m, con, obj, opts);
    end
    
    % Evaluate G at the new parameter sets
    Gs = zeros(nPair, nSam, nSam);
    for iPair = 1:nPair
        if verbose; fprintf('Mapping parameters (#%i) and (#%i)...\n', pairIndexes(iPair,1), pairIndexes(iPair,2)); end
        if opts.PlotProgress;
            % Prepare figure
            iFig = figure();
        end
        
        % Loop through all enumerated values
        for i = 1:nSam
            for j = 1:nSam
                % Fetch the parameters for this point
                knew = paramSets{iPair, i, j}(1:m.nk);
                x0new = paramSets{iPair, i, j}(m.nk+(1:m.nx));
                
                % Evaluate enumerated value
                if strcmpi(opts.EnumSource, 'FIM') || strcmpi(opts.EnumSource, 'hessian')
                    T0 = [m.k(opts.UseParams); m.x0(opts.UseICs)];
                    T = [knew(opts.UseParams); x0new(opts.UseICs)];
                    if opts.Normalized
                        % Use log parameters
                        deltaT = log(T0 ./ T);
                    else
                        % Use regular parameters
                        deltaT = T0 - T;
                    end
                    
                    % Use Taylor expansion to estimate
                    if strcmp(opts.EnumSource, 'FIM')
                        Gs(iPair, i, j) = G0 + deltaT.' * H * deltaT;
                    else
                        Gs(iPair, i, j) = G0 + 0.5 * deltaT.' * H * deltaT;
                    end
                else
                    % Update model to reflect this points parameters
                    mnew = m.Update(knew, x0new, m.q);
                    
                    % Get the objective function value at this point
                    Gs(iPair, i, j) = ObjectiveValue(mnew, con, obj, opts);
                end
                
                % Plot progress
                if opts.PlotProgress
                    figure(iFig)
                    hold off
                    % Pull out the current slice and subtract the min to normalize to 0
                    scaledG(1:nSam,1:nSam) = Gs(iPair,:,:) - min(min(Gs(iPair,:,:)));
                    
                    % Compute the percentile cutoffs
                    percentSpan = linspace(0, 100, nSam^2).';
                    percentile = prctile(vec(scaledG), percentSpan).';
                    
                    % Compute position of each goal along the percentile
                    scaledG = piecewiselinear(scaledG, percentile, percentSpan);
                    scaledG = reshape(scaledG, nSam, nSam);
                    
                    % Draw matrix as image
                    if strcmp(opts.PlotType, 'heat')
                        imagesc(flipud(scaledG.')) % Images have origin in upper left with vertical x and horizontal y
                    else% opts.PlotType == 'contour'
                        contour(scaledG.') % Contours have origin in lower left with vertical x and horizontal y
                    end
                    set(gca, 'Position', [0,0,1,1])
                    set(gca, 'XTick', [])
                    set(gca, 'YTick', [])
                end
            end
        end
    end
    
    % Plot results
    if opts.ToPlot || nargout == 0
        figure()
        for iPair = 1:nPair
            % Pull out the current slice and subtract the min to normalize to 0
            scaledG(1:nSam,1:nSam) = Gs(iPair,:,:) - min(min(Gs(iPair,:,:)));
            
            % Compute the percentile cutoffs
            percentSpan = linspace(0, 100, nSam^2).';
            percentile = prctile(vec(scaledG), percentSpan).';
            
            % Compute position of each goal along the percentile
            scaledG = piecewiselinear(scaledG, percentile, percentSpan);
            scaledG = reshape(scaledG, nSam, nSam);
            
            % Draw matrix as image
            axes()
            if strcmp(opts.PlotType, 'heat')
                imagesc(flipud(scaledG.')) % Images have origin in upper left with vertical x and horizontal y
            else% opts.PlotType == 'contour'
                contour(scaledG.') % Contours have origin in lower left with vertical x and horizontal y
            end
            width = 1 / (nT - 1);
            set(gca, 'ActivePositionProperty', 'outerposition')
            set(gca, 'OuterPosition', [(pairIndexes(iPair,2)-2)*(width), 1-pairIndexes(iPair,1)*(width), width, width])
            set(gca, 'Position', [(pairIndexes(iPair,2)-2)*(width), 1-pairIndexes(iPair,1)*(width), width, width])
            set(gca, 'XTick', [])
            set(gca, 'YTick', [])
        end
    end

    % Work-down
    if nargout >= 1
        % Independent Variables
        data.Pairs          = pairIndexes;
        data.ParameterSets  = paramSets;
        % Results
        data.GoalValues     = Gs;
        
        varargout{1} = data;
    end
end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Optimization functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function val = optFun(l)
        % uses global variables: m, con, obj, nTk, nTk, nT, G0, goalCutoffs, opts, opts
        % also inherits variables v and iCut
        
        % Set new parameters
        newK = m.k;
        newIC = m.x0;
        newP = [m.k(opts.UseParams); m.x0(opts.UseICs)];
        if opts.Normalized
            newP = newP + (newP .* v) * l; % v is in log p space
        else
            newP = newP + v * l;
        end
        newK(opts.UseParams) = newP(1:nTk);
        newIC(opts.UseICs) = newP(nTk + (1:nTx));
        curm = m.Update(newK, newIC, m.q);
        
        % Current objective function
        curG = ObjectiveValue(curm, con, obj, opts);
        
        % Difference from baseline
        deltaG = curG - G0;
        
        % Difference from the cutoff we want to reach
        val = deltaG - goalCutoffs(iCut);
    end

%     function curl = walk()
%         % walk until G target is met or p goes negative
%         deltal = d / 10; % step size
%         G = -inf;
%         step = 0;
%         newP = [m.k(opts.UseParams); m.x0(opts.UseICs)];
%         while (G < 0 && all(newP > 0))
%             step = step + 1;
%             
%             % Initialize
%             newK = m.k;
%             newIC = m.x0;
%             newP = [m.k(opts.UseParams); m.x0(opts.UseICs)];
%             
%             % Compute length
%             curl = step * deltal;
%             
%             % Add step
%             newP = newP + v * curl;
%             newK(opts.UseParams) = newP(1:nTk);
%             newIC(opts.UseICs) = newP(nTk + 1:nTx);
%             curm = m.update(newK, newIC);
%             
%             if any(newP <= 0); break; end
%             
%             G = ObjectiveValue(curm, con, obj, opts) - G0 - goalCutoffs(iCut);
%         end
%     end

end