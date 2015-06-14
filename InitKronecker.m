function InitKronecker()
% This script adds all the relevant Kronecker folders to the current path.
% Run this script before trying to use Kronecker or you will be sad.

kroneckerPath = fileparts(mfilename('fullpath'));

% Universal paths
addpath([kroneckerPath '/Source']);
addpath([kroneckerPath '/External']);
addpath([kroneckerPath '/External/ode15sf']);
addpath([kroneckerPath '/External/libSBML-5.11.4-matlab']);

% Compatibility paths
% Add our version if Matlab version does not exist
compatibility_files = {'assert', 'bsxfun', 'ismatrix', 'padarray', 'strjoin', 'histcounts'};

for i = 1:numel(compatibility_files)
    if ~(exist(compatibility_files{i}, 'builtin') || exist(compatibility_files{i}, 'file'))
        addpath([kroneckerPath '/Compatibility/' compatibility_files{i}])
    end
end

disp('KroneckerBio v0.5.0');