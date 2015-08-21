function tf = issym(obj)
%issym True for symbolic arrays.
%   issym returns true if obj is a symbolic array and false otherwise.

% (c) 2015 David R Hagen
% This work is released under the MIT license.

tf = isa(obj, 'sym');
