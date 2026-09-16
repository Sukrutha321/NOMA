function genBinomCoeff = nchoosek_sVer2(n,s)
% clc
% clear 
% n = 3
% s = 8
prime_vector = ones(1,s);
if (n == 0)
    genBinomCoeff = 1;
elseif (n == 1)
    genBinomCoeff = prime_vector;
elseif (n > 1)
    recurs_vector = prime_vector;
    for i = 2:n
        recurs_vector = conv(recurs_vector,prime_vector);
    end
    genBinomCoeff = recurs_vector;
end
% genBinomCoeff