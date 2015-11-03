diff_vectorized_simplifyFraction := proc(nums, dens)
begin
    if testtype(nums,Dom::Matrix()) then
        [nnonzeroterms,tempnottobeused] := linalg::matdim(nums);
    else
        nnonzeroterms := 1;
    end_if;
    DerivativeVector := matrix(nnonzeroterms, 1);
    if nnonzeroterms = 1 then
        DerivativeVector[1,1] := simplifyFraction(diff(nums,dens));
    else
        for nonzerotermindex from 1 to nnonzeroterms do
            DerivativeVector[nonzerotermindex,1] := simplifyFraction(diff(nums[nonzerotermindex,1],dens[nonzerotermindex,1]));
        end_for;
    end_if;
    return(DerivativeVector)
end_proc