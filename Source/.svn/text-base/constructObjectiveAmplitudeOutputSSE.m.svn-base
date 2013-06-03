function obj = constructObjectiveAmplitudeOutputSSE(m, outputIndex, data, eventIndex, tOn)
% 
% 
% function obj = constructObjectiveOutputAmplitudeOutputSSE(m, outputIndex, data)
% 
% This function constructs an objective function that returns the value of
% the timing of an event, as well as its derivatives.
% 
% Inputs
%       m               -   the model from which output is collected.
%       outputIndex     -   the index of the output whose peak was
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
%   constructObjectiveTimingSpeciesSSE
%   constructObjectiveTimingOutputValue

% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

assert(length(outputIndex) == 1, '', 'outputIndex must be of length 1.');

obj.update      = @update;
obj.G           = @G;
obj.dGdx        = @dGdx;
obj.dGdp        = @dGdp;
obj.d2Gdx2      = @d2Gdx2;
obj.d2Gdp2      = @d2Gdp2;
obj.d2Gdxdp     = @d2Gdxdp;
obj.d2Gdpdx     = @d2Gdpdx;

%% update function
    function objNew = update(mNew)
        objNew = constructObjectiveAmplitudeOutputSSE(mNew, outputIndex, data, eventIndex);
    end

%% G
    function [val stopTimes] = G(t, xSol, u, nSim)
        if isempty(data.i), val = 0; stopTimes = []; return, end
        % (1) Event detection
        val = inf;
        stopTimes = [];
        if isempty(xSol(nSim).ie), return, end      % No events detected
        xSol(nSim) = processp53Peaks(m, xSol(nSim), tOn, data.order);

        % ***************************************************
        % Based on data.i containing the indices of the pulses that are interesting
        % ***************************************************
        ei = find(xSol(nSim).ie == eventIndex);
        if length(ei) < length(data.i), return, end    % Not enough events detected
        
        % (2) Computing objective function value
        val       = (m.c(outputIndex, :)*xSol(nSim).ye(1:m.nX, ei(1:length(data.i))) - data.y);
        val       =  val*val';
        stopTimes =  xSol(nSim).xe(ei(1:length(data.i)));
    end

%% dGdx
    function val = dGdx(t, xSol, u, nSim)
        if isempty(data.i), val = zeros(1, m.nX); return, end
        val = zeros(1, m.nX);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        xSol(nSim) = processp53Peaks(m, xSol(nSim), tOn, data.order);

        % ***************************************************
        % Based on data.i containing the indices of the pulses that are interesting
        % ***************************************************
        ei = find(xSol(nSim).ie == eventIndex);
        if length(ei) < length(data.i), return, end    % Not enough events detected
        
        currentTime = find(t == xSol(nSim).xe(ei(1:length(data.i))));
        if ~isempty(currentTime)
            xe      = deval(xSol(nSim), t, 1:m.nX);
            val       = 2*(m.c(outputIndex, :)*xe - data.y(:,currentTime)).'*m.c(outputIndex,:);
        end
    end

%% dGdp
    function val = dGdp(t, xSol, u, nSim)
        val = zeros(1, m.nP);
    end

%% d2Gdx2
    function val = d2Gdx2(t, xSol, u, nSim)
        val = zeros(m.nX, m.nX);
    end

%% d2Gdp2
    function val = d2Gdp2(t, xSol, u, nSim)
        val = zeros(m.nP, m.nP);
    end

%% d2Gdxdp
    function val = d2Gdpdx(t, xSol, u, nSim)
        val = zeros(m.nX, m.nP);
    end

%% d2Gdpdx
    function val = d2Gdxdp(t, xSol, u, nSim)
        val = zeros(m.nP, m.nX);
    end
end