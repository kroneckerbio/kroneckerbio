fastchar := proc(f)
begin
    //return(coerce(f,DOM_STRING))
    numberofterms := nops(f);
    stringtoreturn := "{";
    if numberofterms = 1 then
        stringtoreturn := stringtoreturn."'".expr2text(f)."'";
    else
        for termindex from 1 to numberofterms do
            if termindex = numberofterms then
                stringtoreturn := stringtoreturn."'".f[termindex,1]."'";
            else
                stringtoreturn := stringtoreturn."'".f[termindex,1]."';";
            end_if;
        end_for;
    end_if;
    stringtoreturn := stringtoreturn."}";
    return(stringtoreturn)
end_proc