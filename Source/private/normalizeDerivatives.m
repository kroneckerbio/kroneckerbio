function [dadT, d2adT2] = normalizeDerivatives(T, dadT, d2adT2)
% [dadT, d2adT2] = normalizeDerivatives(T, dadT, d2adT2)
% Convert sensitivities da/dT and/or curvatures d2a/dT2 at various times t
% into da/d(logT) or d2a/d(logT)2. a can be x (the states), u (the inputs),
% or y (the outputs). T is the list of used parameters.
%
% Input arguments:
%   T:      Vector of parameters chosen by UseParams, UseSeeds,
%           UseInputControls, and UseSeedControls
%   dadT:   na*nT-by-nt matrix of sensitivities of a with respect to T. na
%           is the number of values in a, nT is the number of used parameters,
%           and nt is the number of time points for which values are
%           provided.
%   d2adT2: (optional) na*nT*nT-by-nt matrix of unnormalized curvatures of a with
%           respect to T.
%
% Output arguments:
%   dadT:   Normalized sensitivities da/dlogT.
%   d2adT2: Normalized curvatures d2a/dlogT2. Only returned if unnormalized
%           curvatures are provided.


order = nargin - 1;

nT = numel(T);
nt = size(dadT, 2);
na = size(dadT,1)/nT;

if order >= 1
    
    T_stack_a = vec(repmat(row(T), na, 1));
    dadT = bsxfun(@times, dadT, T_stack_a);
    
    if order >= 2
        
        assert(size(d2adT2,1) == nT.^2.*na, 'KroneckerBio:normalizeDerivatives:d2adT2RowSizeMismatch', 'd2adT2 did not have na*nT*nT rows. This is probably a Kronecker bug.')
        assert(size(d2adT2,2) == nt, 'KroneckerBio:normalizeDerivatives:d2adT2ColumnSizeMismatch', 'd2adT2 did not have nt columns. This is probably a Kronecker bug.')
        
        TT_stack_a = vec(repmat(row(T), na*nT,1)) .* vec(repmat(row(T), na,nT));
        d2adT2 = bsxfun(@times, d2adT2, TT_stack_a) + getNormalizationTerm(dadT);
        
    end
    
end

    function dadT_diag_allt = getNormalizationTerm(dadT)
        % This function performs diagonal expansion for each time point's
        % sensitivities, then turns the result into a column vector, which
        % is stored in the respective time column of an (na*nT*nT)-by-nt
        % matrix. The resulting matrix matches the size of the curvature
        % matrices.
        dadT_diag_allt = zeros(na*nT*nT, nt);
        for ti = 1:nt
            dadT_diag = diagonalExpansion(dadT(:,ti));
            dadT_diag_allt(:,ti) = dadT_diag(:);
        end
    end


    function dadT_diag = diagonalExpansion(dadT)
        % "Diagonally expanding" an na-by-nT matrix dadT means making an
        % (na*nT)-by-nT matrix like so:
        %
        %   [dadT(:,1) 0 0 .  . . 0
        %    0 dadT(:,2) 0
        %    0  0 dadT(:,3)
        %    .              .
        %    .                .
        %    0                  dadT(:,nT)]
        %
        % where 0 indicates an na-by-1 vector of zeros.
        
        dadT_linind = sub2ind([na nT nT], repmat((1:na).',nT,1), vec(repmat(1:nT, na, 1)), vec(repmat(1:nT, na, 1)));
        dadT_diag = zeros(na*nT, nT);
        dadT_diag(dadT_linind) = dadT(:);
        
    end



end