function InitKronecker()
% This script adds all the relevant Kronecker folders to the current path.
% Run this script before trying to use Kronecker or you will be sad.

kroneckerPath = fileparts(mfilename('fullpath'));

% Universal paths
addpath([kroneckerPath '/Source']);
addpath([kroneckerPath '/External']);
addpath([kroneckerPath '/External/ode15sf']);

isoctave = exist('OCTAVE_VERSION', 'builtin');

if ~isoctave
    % Compatibility paths
    matlabVersion = version('-release');
    year = str2double(regexp(matlabVersion, '\d*', 'match', 'once'));
    
    if year <= 2006
        addpath([kroneckerPath '/Compatibility/2006'])
    end
else%isoctave
end

disp('KroneckerBio v0.3.0 alpha');