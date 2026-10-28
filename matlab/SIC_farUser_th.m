function BERth = SIC_farUser_th(SNR_dB,a,M,T,channel,m_factor)
% clc
% clear
% PL = 10.^([0 6]./10);
% SNR_dB = 0:35;
% M = [256 256];
% a = equal_distance(log2(M)).';
% T = [1 2^8];
if (sum(M==[2,2])==2)
    M = [4,4];
    kappa = [1 1];
else
    kappa = 2/3*(M-1);
end
maxT = max(T);
minT = min(T);
tau = maxT./T;
lambda = sqrt(M);
beta = sqrt(M(1)^tau(1)*M(2))/2;
if (beta == inf)
    beta_log = tau(1)/2*log10(M(1)) + 1/2*log10(M(2)) - log10(2);
    beta = 10.^(beta_log);
end

% genBinomCoeff = nchoosek_s(T(2),sqrt(M(1)));
genBinomCoeff = nchoosek_sVer2(T(2),sqrt(M(1)));
% BERth = 0;
for k = 1:1:log2(sqrt(M(2)))
    BERth_k{k} = 0;
    for m = 1:1:(1 - 2^-k)*sqrt(M(2))
        for i = 0:1:tau(1)*(sqrt(M(1)) - 1)
            c1{k}(m,i+1) = genBinomCoeff(i+1);
            c2{k}(m,i+1) = 2^(k-1) - floor((m-1)*2^(k-1)/sqrt(M(2)) + 1/2);
            c3{k}(m,i+1) = (-1)^floor((m-1)*2^(k-1)/sqrt(M(2)));
            delta_i{k}(m,i+1) = (2*i - (sqrt(M(1)) - 1)*tau(1));
            delta_m{k}(m,i+1) = (2*m-1)*tau(1);
%             BERth_k{k} = BERth_k{k} + c1*c2*c3*qfunc(A*sqrt(gamma/T(2)));
        end
    end
    c_k{k} = c1{k}.*c2{k}.*c3{k};
%     delta_i{k}
%     delta_m{k}
    A{k} = delta_i{k}.*sqrt(a(1)/kappa(1)) + delta_m{k}.*sqrt(a(2)/kappa(2));
%     BERth_k{k} = BERth_k{k}./beta;
%     BERth = BERth + BERth_k{k};
end
%%
% BERth = 2/log2(M(2)).*BERth
if(strcmp(channel,'AWGN'))
    PL = 10.^([0]./10);
else
    PL = 10.^([6]./10);
end
SNR = 10.^(SNR_dB./10);
Sd = sqrt(PL)./sqrt(2.*SNR);
BERth = 0;
for k = 1:1:log2(sqrt(M(2)))
    delta = A{k};
    c{k} = c_k{k}(:);
    gamma{k} = (delta(:)./Sd).^2./maxT;
    if(strcmp(channel,'Ray'))
        PHI{k} = 0.5.*(1 - sqrt(gamma{k}./(2 + gamma{k})));
        BERth_k{k} = 1./beta.*sum(c{k}.*PHI{k},1);
    elseif(strcmp(channel,'Nak'))
        mu{k} = sqrt(gamma{k}./(2.*m_factor(2) + gamma{k}));
        k_nak = [0:1:m_factor(2)-1].';
        k_nak = reshape(k_nak,[1,1,length(k_nak)]);
        PHI{k} = mu{k}.*sum(round(exp(gammaln(2.*k_nak+1)-gammaln(k_nak+1)-gammaln(k_nak+1))).*((1 - mu{k}.^2)./4).^k_nak,3);
        BERth_k{k} = 1./beta.*sum(c{k}.*0.5.*(1 - PHI{k}),1);
    elseif(strcmp(channel,'Rice'))
        K_dB = m_factor(2);
        K_r = 10.^(K_dB./10);
        gamma{k} = (delta(:)./Sd).^2./maxT/(1+K_r);
        a_rice{k} = sqrt(K_r./(2.*(1+0.5.*gamma{k}))).*sqrt(1+2.*0.5.*gamma{k}-2.*sqrt(0.5.*gamma{k}.*(1+0.5.*gamma{k})));
        b_rice{k} = sqrt(K_r./(2.*(1+0.5.*gamma{k}))).*sqrt(1+2.*0.5.*gamma{k}+2.*sqrt(0.5.*gamma{k}.*(1+0.5.*gamma{k})));
        BERth_k{k} = 1./beta.*sum(c{k}.*(marcumq(a_rice{k},b_rice{k}) - 0.5*(1 + sqrt(0.5*gamma{k}./(1 + 0.5*gamma{k}))).*exp(-0.5*(a_rice{k}.^2 + b_rice{k}.^2)).*besseli(0,a_rice{k}.*b_rice{k})));
    elseif(strcmp(channel,'AWGN'))
%         BERth_k{k} = 1./beta*sum(c{k}.*qfunc(sqrt(gamma{k})),1);
        BERth_k{k} = 1./2./beta*sum(c{k}.*erfc(sqrt(gamma{k}./2)),1);
    end
    BERth = BERth + BERth_k{k};
end
BERth = 2/log2(M(2)).*BERth;
