function fstr = fastchar(f)
% Convert symbolic array to cell array of strings
%
% For older versions of Matlab, this was a fast conversion function, but now it
% is a simple wrapper.

fstr = cellstr(string(f));
