diff_vectorized := proc(nums, dens)
begin
    [nzterms,temp] := linalg::matdim(nums);
    Der := matrix(nzterms, 1);
    for nti from 1 to nzterms do
        Der[nti,1] := simplifyFraction(diff(nums[nti,1],dens[nti,1]));
    end_for;
    return(Der)
end_proc