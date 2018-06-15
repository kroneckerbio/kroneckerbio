function InitKronecker()
% This script adds all the relevant Kronecker folders to the current path.
% Run this script before trying to use Kronecker or you will be sad.

kroneckerPath = fileparts(mfilename('fullpath'));

% Universal paths
addpath([kroneckerPath '/Source']);
addpath([kroneckerPath '/External']);
addpath([kroneckerPath '/External/ode15sf']);
addpath([kroneckerPath '/External/libSBML-5.17.0-matlab']);

% Compatibility paths
% Add our version if Matlab version does not exist
compatibility_files = {'strjoin'};

for i = 1:numel(compatibility_files)
    if ~(exist(compatibility_files{i}, 'builtin') || exist(compatibility_files{i}, 'file'))
        addpath([kroneckerPath '/Compatibility/' compatibility_files{i}])
    end
end

% Grab and display version
VERSION = KroneckerVersion();
[status, cmdout] = system('git rev-parse HEAD');
if status == 0
    COMMIT = [' based on commit ' cmdout(1:7)];
else
    COMMIT = '';
end
disp(['KroneckerBio v' VERSION COMMIT]);
