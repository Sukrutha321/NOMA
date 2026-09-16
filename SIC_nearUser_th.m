function BERth = SIC_nearUser_th(SNR_dB,a,M,T,channel,m_factor)
% clc
% clear
% close all
% T = [1 2];
M = [2 2];
% a = [0.2 0.8];
% SNR_dB = 0:50;
SNR = 10.^(SNR_dB./10);
% m_factor = [2 2];
maxT = max(T);
minT = min(T);
keyA1 = [-T(2):2:T(2)];
keyA2 = [T(2).*ones(1,length(keyA1))];
keyA = keyA1.*sqrt(a(1)) + keyA2.*sqrt(a(2));
A10 = sqrt(a(1));
A12 = [1 2]*sqrt(a).';
A_12 = [-1 2]*sqrt(a).';
A_row1 = A10;
A_row2 = keyA(1:end-1)./sqrt(T(2));
A_row3 = keyA(2:end)./sqrt(T(2));
A_row4 = keyA(1:end-1)./sqrt(T(2));
A_row5 = keyA./sqrt(T(2));
coeffKey1 = nchoosek_sVer2(T(2)-1,2);
coeffKey2 = nchoosek_sVer2(T(2),2);
c_row1 = 1;
c_row2 = 1./2^(T(2)).*coeffKey1;
c_row3 = 1./2^(T(2)).*coeffKey1;
c_row4 = -1./2^(T(2)).*coeffKey1;
c_row5 = -1./2^(T(2)).*coeffKey2;

% BER_row1 = c_row1.*qfunc(sqrt(gamma.*A_row1.^2));
% BER_row2 = sum([c_row2].'.*qfunc(sqrt(gamma.*[A_row2.^2].')),1);
% BER_row3 = sum([c_row3].'.*qfunc(sqrt(gamma.*A12.^2)).*qfunc(sqrt(gamma.*[A_row3.^2].')),1);
% BER_row4 = sum([c_row4].'.*qfunc(sqrt(gamma.*A_12.^2)).*qfunc(sqrt(gamma.*[A_row4.^2].')),1);
% BER_row5 = sum([c_row5].'.*qfunc(sqrt(gamma.*A10.^2)).*qfunc(sqrt(gamma.*[A_row5.^2].')),1);

BER_row1 = 0;
BER_row2 = 0;
BER_row3 = 0;
BER_row4 = 0;

PL = 10.^([0]./10);

SNR = 10.^(SNR_dB./10);
Sd = sqrt(PL)./sqrt(2.*SNR);

if(strcmp(channel,'Nak')|strcmp(channel,'Ray'))
    if(strcmp(channel,'Ray'))
        m_factor = [1 1];
    end
    BER_row1 = c_row1.*Epsi(A_row1,SNR_dB,m_factor(1),PL(1));
    for i = 1:length(c_row2)
        BER_row2 = BER_row2 + c_row2(i).*Epsi(A_row2(i),SNR_dB,m_factor(1),PL(1));
        BER_row3 = BER_row3 + c_row3(i).*Phy(A12,A_row3(i),SNR_dB,m_factor(1),PL(1));
        BER_row4 = BER_row4 + c_row4(i).*Phy(A_12,A_row4(i),SNR_dB,m_factor(1),PL(1));
    end
    
    BER_row5 = 0;
    for i = 1:length(c_row5)
        BER_row5 = BER_row5 + c_row5(i).*Phy(A10,A_row5(i),SNR_dB,m_factor(1),PL(1));
    end
elseif(strcmp(channel,'AWGN'))
    gamma = (1./Sd).^2;
    BER_row1 = c_row1.*qfunc(A_row1*sqrt(gamma));
    for i = 1:length(c_row2)
        BER_row2 = BER_row2 + c_row2(i).*qfunc(A_row2(i)*sqrt(gamma));
        BER_row3 = BER_row3 + c_row3(i).*qfunc(A12*sqrt(gamma)).*qfunc(A_row3(i)*sqrt(gamma./maxT));
        BER_row4 = BER_row4 + c_row4(i).*qfunc(A_12*sqrt(gamma)).*qfunc(A_row4(i)*sqrt(gamma./maxT));
    end
    
    BER_row5 = 0;
    for i = 1:length(c_row5)
        BER_row5 = BER_row5 + c_row5(i).*qfunc(A10*sqrt(gamma)).*qfunc(A_row5(i)*sqrt(gamma./maxT));
    end
end

    BERth = BER_row1 + BER_row2 + BER_row3 + BER_row4 + BER_row5;
