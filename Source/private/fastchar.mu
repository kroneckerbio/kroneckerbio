fastchar := proc(f)
begin
    //return(coerce(f,DOM_STRING))
    nf := nops(f);
    fstr := "{";
    if nf = 1 then
        fstr := fstr."'".expr2text(f)."'";
    else
        for fi from 1 to nf do
            if fi = nf then
                fstr := fstr."'".f[fi,1]."'";
            else
                fstr := fstr."'".f[fi,1]."';";
            end_if;
        end_for;
    end_if;
    fstr := fstr."}";
    return(fstr)
end_proc