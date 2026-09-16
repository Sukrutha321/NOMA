% JMLSD and SIC comparison
clc
clear
close all
SNR_dB = 0:5:35;
M = [2 2];
a1max = alpha1Max2UE(log2(M));
a = equal_distance(log2(M)).'
T = [1 2];
Iter = 1e5;
channel = 'AWGN';
m_factor = [];
%[BER_JMLSD,~] = JMLSD_sim(SNR_dB,a,M,T,Iter,channel,m_factor);
%BER_SIC = SIC_sim(SNR_dB,a,M,T,Iter,channel,m_factor);
SNR_dB = 0:1:35;
[BERth_JMLSD,~] = JMLSD_th(SNR_dB,a,M,T,channel,m_factor);
BERth_SIC(1,:) = SIC_nearUser_th(SNR_dB,a,M,T,channel,m_factor);
BERth_SIC(2,:) = SIC_farUser_th(SNR_dB,a,M,T,channel,m_factor);
save('/results/results')