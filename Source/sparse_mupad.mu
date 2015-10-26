sparse_mupad := proc(i, j, s, m, n)
begin
    if not testtype(s,Dom::Matrix()) then
        numberofnonzeroentries := 1;
    else
        [numberofnonzeroentries,throwawayvariable] := linalg::matdim(s);
    end_if;
    newmatrixtabletobeusedforinitialization := table();
    if numberofnonzeroentries = 1 then
        newmatrixtabletobeusedforinitialization[i,j] := s;
    else
        for nonzeroindex from 1 to numberofnonzeroentries do
            newmatrixtabletobeusedforinitialization[i[nonzeroindex,1],j[nonzeroindex,1]] := s[nonzeroindex,1];
        end_for;
    end_if;
    //SparseMatrix := matrix(m,n,newmatrixtable);
    return(matrix(m,n,newmatrixtabletobeusedforinitialization))
end_proc