sparse_mupad := proc(i, j, s, m, n)
begin
    [numberofnonzeroentries,throwawayvariable] := linalg::matdim(s);
    newmatrixtabletobeusedforinitialization := table();
    for nonzeroindex from 1 to numberofnonzeroentries do
        newmatrixtabletobeusedforinitialization[i[nonzeroindex,1],j[nonzeroindex,1]] := s[nonzeroindex,1];
    end_for;
    //SparseMatrix := matrix(m,n,newmatrixtable);
    return(matrix(m,n,newmatrixtabletobeusedforinitialization))
end_proc