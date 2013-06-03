function Finv = infoinv(F)

% Desparsify
if issparse(F)
    sparseF = true;
    F = full(F);
else
    sparseF = false;
end

% Size
nT = size(F,1);

% Find specials: zeros and infs
zeroInd = find(diag(F) == 0);
infInd  = find(diag(F) == inf);

% Invert non-special portion of F
invertableInd = true(nT,1);
invertableInd([zeroInd;infInd]) = false;
F = posdef(F(invertableInd,invertableInd));
F = inv(F);

% Create the diagonal of Finv and replace all previous zeros with inf
Finv = zeros(nT,1);
Finv(zeroInd) = inf;
Finv = diag(Finv); % Build the matrix from the diagonal

% Insert all the invertable portions
Finv(invertableInd,invertableInd) = F;

if sparseF
    Finv = sparse(Finv);
end