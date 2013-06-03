function obj = constructObjectiveTimingSpeciesSSE(m, speciesIndex, data)
% 
% 
% function obj = constructObjectiveOutputTimingSpeciesSSE(m, speciesIndex, data)
% 
% This function constructs an objective function that returns the value of
% the timing of an event, as well as its derivatives.
% 
% Inputs
%       m               -   the model from which output is collected.
%       speciesIndex    -   the index of the species whose peak was
%                           detected by the event described in eventIndex.
%       data            -   a structure containing the data with which to
%                           compare.  This structure has two fields:
%                               data.i is a vector of event numbers (i.e. 
%                                      which pulses are of interest)
%                               data.t is the timing of each of these events
% 
% Outputs
%       obj             -   an objective function structure containing
%                           three functions:
%                           obj.G is the objective function value
%                           obj.dGdx is the derivative with respect to species
%                           obj.dGdp is the derivative with respect to parameters
% 
% See also:
%   constructObjectiveOutputValue
%   constructObjectiveOutputISE
%   constructObjectiveOutputSSE
%   constructObjectiveSpeciesValue
%   constructObjectiveSpeciesISE
%   constructObjectiveSpeciesSSE
%   constructObjectiveTimingSpeciesValue
%   constructObjectiveTimingOutputValue
%   constructObjectiveTimingOutputSSE

% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.


assert(length(speciesIndex) == 1, '', 'speciesIndex must be of length 1.');

% Compute the matrix that selects out the interesting species.
C = eye(m.nX);
C = C(speciesIndex,:);

obj.G           = @G;
obj.dGdx        = @dGdx;
obj.dGdp        = @dGdp;

    function [val stopTimes] = G(t, xSol, u, nSim)
        % (1) Event detection
        val = inf;
        stopTimes = [];
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < max(data.i), return, end    % Not enough events detected
        
        % (2) Computing objective function value
        val       = (xSol(nSim).xe(ei(data.i)) - data.t(i));
        val       =  val*val';
        stopTimes =  xSol(nSim).xe(ei(data.i));
    end

    function val = dGdx(t, xSol, u, nSim)
        val = zeros(1, m.nX);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < max(data.i), return, end    % Not enough events detected
        
        currentTime = find(t == xSol(nSim).xe(ei(data.i)));
        if ~isempty(currentTime)
            val     = 2*(t - data.t(currentTime));
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
            dfidx   = C(speciesIndex,:)*dfdx;
            val = val*(-1/(dfidx*f)*dfidx);
        end
    end

    function val = dGdp(t, xSol, u, nSim)
        val = zeros(1, m.nP);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < max(data.i), return, end    % Not enough events detected
        
        currentTime = find(t == xSol(nSim).xe(ei(data.i)));
        if ~isempty(currentTime)
            val     = 2*(t - data.t(currentTime));
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
            dfidx   = C(speciesIndex,:)*dfdx;
            dfdp    = m.dfdp(t, ye, u);
            dfidp   = C(speciesIndex,:)*dfdp;
            val = val*(-1/(dfidx*f)*dfidp);
        end
    end
end