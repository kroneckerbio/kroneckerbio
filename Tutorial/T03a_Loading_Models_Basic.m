%% T03a Basic Loading Models
% Load a model from a kroneckerbio format file
%
% The file format can be found in the comments in Testing.txt
%
% Note that the filename is given as a full path in this example. Only the
%   relative path is actually needed normally, but the full path is needed
%   by the automatic documentation generator.

current_path = fileparts(mfilename('fullpath'));

m = LoadModelMassAction(fullfile(current_path, 'Testing.txt'))
