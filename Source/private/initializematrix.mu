initializematrix := proc(oldmatrix, oldindices, M, N, newindices)
begin
    [nz,n2] := linalg::matdim(oldindices);
    newmatrixtable := table();
    for nzi from 1 to nz do
        newmatrixtable[newindices[nzi,1],newindices[nzi,2]] := oldmatrix[oldindices[nzi,1],oldindices[nzi,2]];
    end_for;
    newmatrix := matrix(M,N,newmatrixtable);
    return(newmatrix)
end_proc: