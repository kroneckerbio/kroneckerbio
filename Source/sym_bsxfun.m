function c = sym_bsxfun(func,a,b)

if isa(a,'sym') || isa(b,'sym')
    as = size(a);
    bs = size(b);
    a_repmat = ones(size(as));
    a_repmat(as == 1) = bs(as == 1);
    b_repmat = ones(size(bs));
    b_repmat(bs == 1) = as(bs == 1);
    
    a = repmat(a, a_repmat);
    b = repmat(b, b_repmat);
    
    c = func(a,b);
else
    c = bsxfun(func,a,b);
end

end