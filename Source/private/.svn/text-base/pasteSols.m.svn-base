function xSol = pasteSols(xSol, xSolNew)
%PASTESOLS Paste in matching values from new solution with defined
%   timepoints

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Paste in new values
for i = 1:length(xSolNew.x)
    tGetInds = find(xSolNew.x(i) == xSol.x);
    xSol.y(:,tGetInds) = repmat(xSolNew.y(:,i),1,length(tGetInds));
end

% If there were events in either, append these as well
if isfield(xSol, 'xe') && isfield(xSolNew, 'xe')
    xSol.xe = [xSol.xe xSolNew.xe];
    xSol.ye = [xSol.ye xSolNew.ye];
    xSol.ie = [xSol.ie xSolNew.ie];
elseif isfield(xSolNew, 'xe')
    xSol.xe = xSolNew.xe;
    xSol.ye = xSolNew.ye;
    xSol.ie = xSolNew.ie;
end