diff_vectorized := proc(nums, dens)
begin
    [nnonzeroterms,tempnottobeused] := linalg::matdim(nums);
    DerivativeVector := matrix(nnonzeroterms, 1);
    for nonzerotermindex from 1 to nnonzeroterms do
        DerivativeVector[nonzerotermindex,1] := simplifyFraction(diff(nums[nonzerotermindex,1],dens[nonzerotermindex,1]));
    end_for;
    return(DerivativeVector)
end_proc