function I = Epsi(x,SNR_dB,m_factor,PL)
    SNR = 10.^(SNR_dB./10);
    Sd = sqrt(PL./(2.*SNR));
    gamma_n = 1./Sd.^2;
    mu = sqrt(gamma_n.*x.^2./(2.*m_factor + gamma_n.*x.^2));
    k = [0:1:m_factor-1].';
    I = 0.5.*(1 - mu.*sum(round(exp(gammaln(2.*k+1)-gammaln(k+1)-gammaln(k+1))).*((1 - mu.^2)./4).^k,1));
