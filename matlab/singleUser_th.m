function BERth = singleUser_th(SNR_dB,a,M,channel,m_factor)
% clc
% clear
% SNR_dB = 0:35;
% 
% M = 4;
% a = 1;
if (M==2)
    M = 4;
    kappa = [1];
else
    kappa = 2/3*(M-1);
end
lambda = sqrt(M);
beta = sqrt(M)/2;
for k = 1:1:log2(sqrt(M))
    BERth_k{k} = 0;
    for m = 1:1:(1 - 2^-k)*sqrt(M)
        c2{k}(m,1) = 2^(k-1) - floor((m-1)*2^(k-1)/sqrt(M) + 1/2);
        c3{k}(m,1) = (-1)^floor((m-1)*2^(k-1)/sqrt(M));
        delta_m{k}(m,1) = (2*m-1);
    end
    c_k{k} = c2{k}.*c3{k};
    A{k} = delta_m{k}.*sqrt(a/kappa);
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
for k = 1:1:log2(sqrt(M))
    delta = A{k};
    c{k} = c_k{k}(:);
    gamma{k} = (delta(:)./Sd).^2;
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
        gamma{k} = (delta(:)./Sd).^2./(1+K_r);
        a_rice{k} = sqrt(K_r./(2.*(1+0.5.*gamma{k}))).*sqrt(1+2.*0.5.*gamma{k}-2.*sqrt(0.5.*gamma{k}.*(1+0.5.*gamma{k})));
        b_rice{k} = sqrt(K_r./(2.*(1+0.5.*gamma{k}))).*sqrt(1+2.*0.5.*gamma{k}+2.*sqrt(0.5.*gamma{k}.*(1+0.5.*gamma{k})));
        BERth_k{k} = 1./beta.*sum(c{k}.*(marcumq(a_rice{k},b_rice{k}) - 0.5*(1 + sqrt(0.5*gamma{k}./(1 + 0.5*gamma{k}))).*exp(-0.5*(a_rice{k}.^2 + b_rice{k}.^2)).*besseli(0,a_rice{k}.*b_rice{k})));
    elseif(strcmp(channel,'AWGN'))
%         BERth_k{k} = 1./beta*sum(c{k}.*qfunc(sqrt(gamma{k})),1);
        BERth_k{k} = 1./2./beta*sum(c{k}.*erfc(sqrt(gamma{k}./2)),1);
    end
    BERth = BERth + BERth_k{k};
end
BERth = 2/log2(M).*BERth;
