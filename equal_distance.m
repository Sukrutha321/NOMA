function a = equal_distance(b)
M = 2.^b; % only rectangular QAM
N = length(M);
kappa = ones(1,N);
Lambda = ones(1,N);
symKey = cell(1,N);
LutKey = cell(1,N);
DataKeyLUT = cell(1,N);
xLength = ones(1,N);
yLength = ones(1,N);

for i = 1:N
    symKey{i} = 0:M(i)-1;
    LutKey{i} = rectGrayQAM_mod(M(i));
    kappa(i) = LutKey{i}*LutKey{i}'/M(i);
    Lambda(i) = 2^(ceil(0.5*log2(M(i))))-1;
    LutKey{i} = LutKey{i}/sqrt(LutKey{i}*LutKey{i}'/M(i));
    xLength(i) = length(unique(real(LutKey{i})));
    yLength(i) = length(unique(imag(LutKey{i})));
end
z = zeros(N,1);
z(N,1) = 1;
Omega = zeros(N,N);
Omega(N,:) = ones(1,N);
for i = 1:N-1
    Omega(i,i) = kappa(i+1)*(Lambda(i)+1)^2;
    Omega(i,i+1) = -kappa(i);
end
a = linsolve(Omega,z);
