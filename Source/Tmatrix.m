function T = Tmatrix(si)

[ro co] = deal(si(1), si(2));

T=sparse(ro*co,co*ro);


ind=[1:ro:ro*co]';

ind=ind-1;
r=1;

for c=1:ro
    ind=ind+1;
    
    for aux=1:length(ind)
    
        T(r,ind(aux))=1;
    r=r+1;
    end

     
end

end
    
     