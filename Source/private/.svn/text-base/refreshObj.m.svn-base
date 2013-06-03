function newObj = refreshObj(m, con, obj, UseParams, UseICs, UseControls)
%REFRESHOBJ Update Kronecker Bio objective functions according to a new
%   model
%
%   Whenever either the topolgy or parameters of a model change, the
%   objective function structures should be updated because they may
%   depend on the changes. This function loops over the objective functions
%   in order to accomplish this common task.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 5
    UseControls = [];
    if nargin < 4
        UseICs = [];
        if nargin < 3
            UseParams = [];
        end
    end
end

nTop = numel(m);
nCon = size(con, 1);
nObj = size(obj, 1);

% Clean-up empties
if isempty(UseParams)
    UseParams = cell(nTop,1);
end
if isempty(UseICs)
    UseICs = cell(nTop,1);
end
if isempty(UseControls)
    UseControls = cell(nCon,nTop);
end

% Make options consistent as cell arrays
if ~iscell(UseParams)
    UseParams = {UseParams};
end
if ~iscell(UseICs)
    UseICs = {UseICs};
end
if ~iscell(UseControls)
    UseControls = {{UseControls}};
end
if nTop > 0 && nCon > 0 && ~iscell(UseControls{1})
    UseControls = {UseControls};
end

newObj = Gzero([nObj,nCon,nTop]);
for iTop = 1:nTop
    for iCon = 1:nCon
        for iObj = 1:nObj
            if size(UseICs{iTop},2) == 1 && size(UseControls,1) == 1 && nCon > 1
                % UseModelICs == true && UseModelInputs == true
                useICsi = UseICs{iTop};
                useControlsi = UseControls{iTop}{iCon};
            elseif size(UseICs{iTop},2) == 1 && nCon > 1
                % UseModelsICs == true && UseModelInputs == false
                useICsi = UseICs{iTop};
                useControlsi = UseControls{iTop}{iCon};
            elseif size(UseControls,1) == 1 && nCon > 1
                % UseModelsICs == false && UseModelInputs == true
                useICsi = UseICs{iTop}(:,iCon);
                useControlsi = UseControls{iTop}{iCon};
            else% (size(UseICs{iTop},2) ~= 1 && size(UseControls,1) ~= 1) || nCon == 1
                % UseModelsICs == false && UseModelInputs == false
                useICsi = UseICs{iTop}(:,iCon);
                useControlsi = UseControls{iTop}{iCon};
            end
            
            % Refresh obj
            newObj(iObj,iCon,iTop) = pastestruct(Gzero(m(iTop)), obj(iObj,iCon,iTop).Update(m(iTop), con(iCon,iTop), UseParams{iTop}, useICsi, useControlsi));
        end
    end
end
