function eve = Event(direction, e)
%Event Contructor for KroneckerBio event struct
%
%   eve = Event(direction, e)
%
%   Inputs
%   direction: [ -1 0 1 ]
%       1 means detect when the function becomes positive, -1 means detect when
%       the function becomes negative, 0 means detect when the function changes
%       sign.
%   e: [ handle @(t,y) ]
%       An function handle of time and output vector ny, returning a numeric 
%       scalar. This is the function which the simulator will try to find zeros
%       of.

% (c) 2015 David R Hagen and Bruce Tidor
% This work is released under the MIT license.


eve.direction = direction;
eve.e = e;
