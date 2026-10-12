function [LutKey,kappa,Lambda] = rectGrayQAM_mod(M)

Lambda = 2^(ceil(0.5*log2(M))) - 1;
moh = 2^(floor(0.5*log2(M))) - 1;

InPhase = -Lambda:2:Lambda;
InPhase = ones(moh+1,Lambda+1) * diag(InPhase);

QuadPhase = [moh:-2:-moh].';

LUT = InPhase + 1j*QuadPhase;

% Binary-to-Gray conversion
symbol = 0:M-1;
symbolGray = bitxor(symbol,floor(symbol/2));

LutKey = LUT(:).';

kappa = LutKey * LutKey' / M;

LutKey = LutKey(symbolGray+1);

end