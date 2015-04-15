sparse_mupad := proc(i, j, s, m, n)
begin
    [nzmax,temp] := linalg::matdim(s);
    newmatrixtable := table();
    for nzi from 1 to nzmax do
        newmatrixtable[i[nzi,1],j[nzi,1]] := s[nzi,1];
    end_for;
    S := matrix(m,n,newmatrixtable);
    return(S)
end_proc