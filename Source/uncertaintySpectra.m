function varargout = uncertaintySpectra(F, labels)
% F is a Fisher information matrix or cell array of them; labels is the
% x-axis labels for each information matarix

% Clean-up inputs
if nargin < 2
    labels = [];
end

% Ensure F is a cell array
if ~iscell(F)
    wasMatrix = true;
    F = {F};
else
    wasMatrix = false;
end

n = numel(F);
nP = length(F{1});

% Compute uncertainty spectra
sigmas = zeros(nP,n);
for i = 1:n
    lambda = infoeig(F{i});
    lambda(lambda < eps) = eps;
    sigmas(:,i) = sqrt(lambda).^-1;
end

% Plot sigmas if nargout is 0
if nargout == 0
    % Find max and min for y-axis
    yMin = 10.^(floor(log10(min(sigmas(:))) - 0.09));
    yMax = 10.^(ceil(log10(max(sigmas(:))) + 0.09));
    
    % Setup figure
    set(gca, 'xtick', 1:n,...
             'xticklabel', labels,...
             'YScale', 'log');
    box on
    xlim([0 n+1])
    ylim([yMin yMax])
    xlabel('Experiment')
    ylabel('Uncertainty')
    
    % Draw lines
    width = 0.5; % Length of lines relative to distance between spectra
    for i = 1:n
        for j = 1:nP
            yVal = sigmas(j,i);
            line([i-width/2, i+width/2], [yVal, yVal], 'Color', 'r');
            
        end
    end
else % Sigmas are wanted, not a plot of them
    if wasMatrix
        varargout{1} = sigmas{1};
    else
        varargout{1} = sigmas;
    end
end
