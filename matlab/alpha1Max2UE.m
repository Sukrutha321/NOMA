function a1max = alpha1Max2UE(b)
M = 2.^b;            %modulation orders vector
N = length(M);       %number of users
kappa = ones(1,N);
Lambda = ones(1,N);
symKey = cell(1,N);
LutKey = cell(1,N);
DataKeyLUT = cell(1,N);
xLength = ones(1,N);
yLength = ones(1,N);

% Generating the constellations for each user
for i = 1:N
    symKey{i} = 0:M(i)-1;
    if (M(i)==4)
        kappa(i) = 2;
        Lambda(i) = 1;
    elseif (M(i)==2)
        kappa(i) = 1;
        Lambda(i) = 1;
    elseif (M(i)==32||128||512||2048||8192||32768)
        LutKey{i} = rectGrayQAM_mod(M(i));
        kappa(i) = LutKey{i}*LutKey{i}'/M(i);
        Lambda(i) = 2^(ceil(0.5*log2(M(i))))-1;
    elseif (M(i)~=4||2)
        LutKey{i} = qammod(symKey{i},M(i));
        kappa(i) = LutKey{i}*LutKey{i}'/M(i);
        Lambda(i) = 2^(ceil(0.5*log2(M(i))))-1;
    end
end

a1max = kappa(1)/(kappa(1) + kappa(2)*Lambda(1)^2);