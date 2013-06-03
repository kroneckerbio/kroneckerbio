function uncertaintyVectors(A, cut, pLabels)

[lambda, Q] = infoeig(A);
nT = numel(lambda);

% Defaults
if nargin < 2
    cut = 0.1;
end
if nargin < 3
    pLabels = cellfun(@num2str, num2cell(1:nT)', 'UniformOutput', false);
end

skip = floor(nT/length(pLabels));

% Sort from largest to smallest
[lambda ind] = sort(lambda, 'Descend');
sigma = lambda.^(-1/2);
Q = Q(:,ind);

imagesc(Q, [-1,1]);
colorbar;
colormap(heatMapColors(256, 1));

axis image;

% Draw gray Gridlines
for r = (1:nT)-0.5
    line([r,  r], [0.5, nT+0.5], 'Color', [1,1,1]*0.2);
    line([0.5, nT+0.5], [r,  r], 'Color', [1,1,1]*0.2);
end

% Draw cutoff lines
evColors = ['g', 'y', 'r'];

for i = 1:numel(cut) 
    ind = find(sigma < cut(i), 1, 'Last');

    if(isempty(ind))
        ind = 0.5;
    else
        ind = ind + 0.5;
    end

    line([ind ind], [0 nT+1], 'Color', evColors(i));
end
    
set(gca, 'xTick', []);
set(gca, 'yTick', 1:skip:nT);
set(gca, 'yTicklabel', pLabels);
set(gca, 'TickLength', [0; 0]);
xlabel('Increasing uncertainty \rightarrow');

end

function map = heatMapColors(N, gamma)

if nargin < 1
    N = 128;
else
    if mod(N,2) == 1
        N = N + 1;
    end
end

if nargin < 2
    gamma = 1;
end

map = [linspace(1, 0, N/2).^gamma,              zeros(1, N/2);
                    zeros(1, N/2), linspace(0, 1, N/2).^gamma;
                                                zeros(1, N)]';
end