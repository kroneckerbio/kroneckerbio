function obj = constructObjectiveTimingOutputSSE(m, outputIndex, data, eventIndex, tOn)
% 
% 
% function obj = constructObjectiveOutputTimingOutputSSE(m, outputIndex, data)
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

nS       = m.nP+m.nX;
nP       = m.nP;
Ix       = eye(m.nX);
Ip       = eye(m.nP);

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
        objNew = constructObjectiveTimingOutputSSE(mNew, outputIndex, data, eventIndex);
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
        % Based on data.i containing the number of pulses that are interesting
        % ***************************************************
%         ei = find(xSol(nSim).ie); % == eventIndex
%         if length(ei) < max(data.i), return, end    % Not enough events detected
%         
%         % (2) Computing objective function value
%         val       = (xSol(nSim).xe(ei(data.i)) - data.t);
%         val       =  val*val';
%         stopTimes =  xSol(nSim).xe(ei(data.i));

        % ***************************************************
        % Based on data.i containing the indices of the pulses that are interesting
        % ***************************************************
        ei = find(xSol(nSim).ie == eventIndex);
        if length(ei) < length(data.i), return, end    % Not enough events detected
        
        % (2) Computing objective function value
        val       = (xSol(nSim).xe(ei(1:length(data.i))) - data.t);
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
        % Based on data.i containing the number of pulses that are interesting
        % ***************************************************
%         ei = find(xSol(nSim).ie); % == eventIndex
%         if length(ei) < max(data.i), return, end    % Not enough events detected
%         
%         currentTime = find(t == xSol(nSim).xe(ei(data.i)));
%         if ~isempty(currentTime)
%             ye      = deval(xSol(nSim), t, 1:m.nX);
%             f       = m.f(t, ye, u);
%             dfdx    = m.dfdx(t, ye, u);
%             dfidx   = m.c(outputIndex,:)*dfdx;
%             dtdx    = -1/(dfidx*f)*dfidx;
%             val     = 2*(t - data.t(currentTime))*dtdx;
%         end
        
        % ***************************************************
        % Based on data.i containing the indices of the pulses that are interesting
        % ***************************************************
        ei = find(xSol(nSim).ie == eventIndex);
        if length(ei) < length(data.i), return, end    % Not enough events detected
        
        currentTime = find(t == xSol(nSim).xe(ei(1:length(data.i))));
        if ~isempty(currentTime)
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
            dfidx   = m.c(outputIndex,:)*dfdx;
            dtdx    = -1/(dfidx*f)*dfidx;
            val     = 2*(t - data.t(currentTime))*dtdx;
        end
    end

%% dGdp
    function val = dGdp(t, xSol, u, nSim)
        if isempty(data.i), val = zeros(1, m.nP); return, end
        val = zeros(1, m.nP);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        xSol(nSim) = processp53Peaks(m, xSol(nSim), tOn, data.order);

        % ***************************************************
        % Based on data.i containing the number of pulses that are interesting
        % ***************************************************
%         ei = find(xSol(nSim).ie); % == eventIndex
%         if length(ei) < max(data.i), return, end    % Not enough events detected
%         
%         currentTime = find(t == xSol(nSim).xe(ei(data.i)));
%         if ~isempty(currentTime)
%             ye      = deval(xSol(nSim), t, 1:m.nX);
%             f       = m.f(t, ye, u);
%             dfdx    = m.dfdx(t, ye, u);
%             dfidx   = m.c(outputIndex,:)*dfdx;
%             dfdp    = m.dfdp(t, ye, u);
%             dfidp   = m.c(outputIndex,:)*dfdp;
%             dtdp    = -1/(dfidx*f)*dfidp;
%             val     = 2*(t - data.t(currentTime))*dtdp;
%         end

        % ***************************************************
        % Based on data.i containing the indices of the pulses that are interesting
        % ***************************************************
        ei = find(xSol(nSim).ie == eventIndex);
        if length(ei) < length(data.i), return, end    % Not enough events detected
        
        currentTime = find(t == xSol(nSim).xe(ei(1:length(data.i))));
        if ~isempty(currentTime)
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
            dfidx   = m.c(outputIndex,:)*dfdx;
            dfdp    = m.dfdp(t, ye, u);
            dfidp   = m.c(outputIndex,:)*dfdp;
            dtdp    = -1/(dfidx*f)*dfidp;
            val     = 2*(t - data.t(currentTime))*dtdp;
        end
    end

%% d2Gdx2
    function val = d2Gdx2(t, xSol, u, nSim)
        val = zeros(m.nX, m.nX);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < max(data.i), return, end    % Not enough events detected
        
        currentTime = find(t == xSol(nSim).xe(ei(data.i)));
        if ~isempty(currentTime)
            ye       = deval(xSol(nSim), t, 1:m.nX);
            f        = m.f(t, ye, u);
            dfdx     = m.dfdx(t, ye, u);
            dfidx    = m.c(outputIndex,:)*dfdx;
            dfdp     = m.dfdp(t, ye, u);
            dfidp    = m.c(outputIndex,:)*dfdp;
            d2fdp2   = m.d2fdp2(t, ye, u);
            d2fidp2  = kron(m.c(outputIndex,:),Ip)*d2fdp2;
            d2fdx2   = m.d2fdx2(t, ye, u);
            d2fidx2  = kron(m.c(outputIndex,:),Ix)*d2fdx2;
            d2fdxdp  = m.d2fdxdp(t, ye, u);
            d2fidxdp = kron(m.c(outputIndex,:),Ix)*d2fdxdp;
            d2fdpdx  = m.d2fdpdx(t, ye, u);
            d2fidpdx = kron(m.c(outputIndex,:),Ip)*d2fdpdx;
            dtdx2    = -1/(dfidx*f)*(d2fidx2  + -1/(dfidx*f)*d2fidx2*f*dfidx  + -1/(dfidx*f)*dfdx'*dfidx'*dfidx) + 1/(dfidx*f)^2*dfidx'*(f'*(-1/(dfidx*f)*d2fidx2*f*dfidx + d2fidx2)  + dfidx*(-1/(dfidx*f)*dfdx*f*dfidx + dfdx));
            dtdx     = -1/(dfidx*f)*dfidx;
            val      = 2*(t - data.t(currentTime))*dtdx2 + 2*dtdx'*dtdx;
        end
    end

%% d2Gdp2
    function val = d2Gdp2(t, xSol, u, nSim)
        val = zeros(m.nP, m.nP);
        
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < max(data.i), return, end    % Not enough events detected
        
        currentTime = find(t == xSol(nSim).xe(ei(data.i)));
        if ~isempty(currentTime)
            ye       = deval(xSol(nSim), t, 1:m.nX);
            f        = m.f(t, ye, u);
            dfdx     = m.dfdx(t, ye, u);
            dfidx    = m.c(outputIndex,:)*dfdx;
            dfdp     = m.dfdp(t, ye, u);
            dfidp    = m.c(outputIndex,:)*dfdp;
            d2fdp2   = m.d2fdp2(t, ye, u);
            d2fidp2  = kron(m.c(outputIndex,:),Ip)*d2fdp2;
            d2fdx2   = m.d2fdx2(t, ye, u);
            d2fidx2  = kron(m.c(outputIndex,:),Ix)*d2fdx2;
            d2fdxdp  = m.d2fdxdp(t, ye, u);
            d2fidxdp = kron(m.c(outputIndex,:),Ix)*d2fdxdp;
            d2fdpdx  = m.d2fdpdx(t, ye, u);
            d2fidpdx = kron(m.c(outputIndex,:),Ip)*d2fdpdx;
            dtdp2    = -1/(dfidx*f)*(d2fidp2  + -1/(dfidx*f)*d2fidpdx*f*dfidp + -1/(dfidx*f)*dfdp'*dfidx'*dfidp) + 1/(dfidx*f)^2*dfidp'*(f'*(-1/(dfidx*f)*d2fidx2*f*dfidp + d2fidxdp) + dfidx*(-1/(dfidx*f)*dfdx*f*dfidp + dfdp));
            dtdp     = -1/(dfidx*f)*dfidp;
            val      = 2*(t - data.t(currentTime))*dtdp2 + 2*dtdp'*dtdp;
        end
    end

%% d2Gdxdp
    function val = d2Gdpdx(t, xSol, u, nSim)
        val = zeros(m.nX, m.nP);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < max(data.i), return, end    % Not enough events detected
        
        currentTime = find(t == xSol(nSim).xe(ei(data.i)));
        if ~isempty(currentTime)
            ye       = deval(xSol(nSim), t, 1:m.nX);
            f        = m.f(t, ye, u);
            dfdx     = m.dfdx(t, ye, u);
            dfidx    = m.c(outputIndex,:)*dfdx;
            dfdp     = m.dfdp(t, ye, u);
            dfidp    = m.c(outputIndex,:)*dfdp;
            d2fdp2   = m.d2fdp2(t, ye, u);
            d2fidp2  = kron(m.c(outputIndex,:),Ip)*d2fdp2;
            d2fdx2   = m.d2fdx2(t, ye, u);
            d2fidx2  = kron(m.c(outputIndex,:),Ix)*d2fdx2;
            d2fdxdp  = m.d2fdxdp(t, ye, u);
            d2fidxdp = kron(m.c(outputIndex,:),Ix)*d2fdxdp;
            d2fdpdx  = m.d2fdpdx(t, ye, u);
            d2fidpdx = kron(m.c(outputIndex,:),Ip)*d2fdpdx;
            d2tdxdp  = -1/(dfidx*f)*(d2fidxdp + -1/(dfidx*f)*d2fidx2*f*dfidp  + -1/(dfidx*f)*dfdx'*dfidx'*dfidp) + 1/(dfidx*f)^2*dfidx'*(f'*(-1/(dfidx*f)*d2fidx2*f*dfidp + d2fidxdp) + dfidx*(-1/(dfidx*f)*dfdx*f*dfidp + dfdp));
            dtdx     = -1/(dfidx*f)*dfidx;
            dtdp     = -1/(dfidx*f)*dfidp;
            val      = 2*(t - data.t(currentTime))*d2tdxdp + 2*dtdx'*dtdp;
        end
    end

%% d2Gdpdx
    function val = d2Gdxdp(t, xSol, u, nSim)
        val = zeros(m.nP, m.nX);
        
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < max(data.i), return, end    % Not enough events detected
        
        currentTime = find(t == xSol(nSim).xe(ei(data.i)));
        if ~isempty(currentTime)
            ye       = deval(xSol(nSim), t, 1:m.nX);
            f        = m.f(t, ye, u);
            dfdx     = m.dfdx(t, ye, u);
            dfidx    = m.c(outputIndex,:)*dfdx;
            dfdp     = m.dfdp(t, ye, u);
            dfidp    = m.c(outputIndex,:)*dfdp;
            d2fdp2   = m.d2fdp2(t, ye, u);
            d2fidp2  = kron(m.c(outputIndex,:),Ip)*d2fdp2;
            d2fdx2   = m.d2fdx2(t, ye, u);
            d2fidx2  = kron(m.c(outputIndex,:),Ix)*d2fdx2;
            d2fdxdp  = m.d2fdxdp(t, ye, u);
            d2fidxdp = kron(m.c(outputIndex,:),Ix)*d2fdxdp;
            d2fdpdx  = m.d2fdpdx(t, ye, u);
            d2fidpdx = kron(m.c(outputIndex,:),Ip)*d2fdpdx;
        	d2tdpdx  = -1/(dfidx*f)*(d2fidpdx + -1/(dfidx*f)*d2fidpdx*f*dfidx + -1/(dfidx*f)*dfdp'*dfidx'*dfidx) + 1/(dfidx*f)^2*dfidp'*(f'*(-1/(dfidx*f)*d2fidx2*f*dfidx + d2fidx2)  + dfidx*(-1/(dfidx*f)*dfdx*f*dfidx + dfdx));
            dtdx     = -1/(dfidx*f)*dfidx;
            dtdp     = -1/(dfidx*f)*dfidp;
            val      = 2*(t - data.t(currentTime))*d2tdpdx + 2*dtdp'*dtdx;
        end
    end
end