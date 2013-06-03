function newCon = refreshCon(m, con)
%REFRESHOBJ Update Kronecker Bio experimental conditions a new model
%
%   Whenever either the topolgy or parameters of a model change, the
%   experimental condition structures should be updated because they may
%   depend on the changes. This function loops over the experimental
%   conditions in order to accomplish this common task.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

nTop = numel(m);
nCon = size(con,1);
newCon = Uzero([nCon,nTop]);
for iTop = 1:nTop
    for iCon = 1:nCon
        newCon(iCon,iTop) = pastestruct(Uzero(m(iTop)), con(iCon,iTop).Update(con(iCon,iTop).x0, con(iCon,iTop).q));
    end
end
