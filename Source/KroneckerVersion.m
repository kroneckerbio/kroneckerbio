function [str, epoch, major, minor, bugfix] = KroneckerVersion()
%KroneckerVersion Return the KroneckerBio version
%
%   [version, epoch, major, minor, bugfix] = KroneckerVersion()
%
%   Outputs
%   str: [ string ] 
%       The version as a human readible text string. (e.g. '1.2.3.4')
%   epoch: [ numeric scalar ]
%       The KroneckerBio epoch number (e.g. 1)
%   major: [ numeric scalar ]
%       The major revision number (e.g. 2)
%   minor: [ numeric scalar ]
%       The minor revision number (e.g. 3)
%   bugfix: [ numeric scalar ]
%       The bugfix or point number (e.g. 4)

% (c) 2016 David Hagen & Joshua Apgar
% This work is released under the MIT license.

kroneckerPath = fileparts(mfilename('fullpath'));

str = strtrim(fileread([kroneckerPath '/../VERSION']));
parts = strsplit(str, '.');
epoch = str2double(parts{1});
major = str2double(parts{2});
minor = str2double(parts{3});
bugfix = str2double(parts{4});
