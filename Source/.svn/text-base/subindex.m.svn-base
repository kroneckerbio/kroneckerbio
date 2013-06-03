function A = subindex(A, varargin)
%SUBINDEX Subscript indexes
%
%   A = subindex(A, varargin)
%
%   Inexplicably, Matlab does not allow subscripting to be chained or used
%   to access the output of a function, e.g.
%
%   pixel = myimage(x1:x2,y1:y2)(index)   % error
%   width = size(A)(2)                    % error & arguably poor example.
%
%   Such functionality can be generally useful for code reability or
%   critically important if an anonymous function handle, which cannot
%   store intermediate variables, needs to be created. This function
%   provides inline access to subscripting.
%
%   Examples:
%   pixel = subindex(myimage(x1:x2,y1:y2), index)
%   width = subindex(size(A), 2)

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

A = A(varargin{:});
end