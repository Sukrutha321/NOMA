clc;
clear;
close all;

%% Parameters
M = [2 2];
T = [1 2];
channel = 'AWGN';
m_factor = [];

%% ============================================================
% FIGURE 6(a): BER versus SNR
% alpha1 = 0.2, alpha2 = 0.8
% ============================================================

SNR_dB = 0:1:15;

a = [0.2 0.8];

% JMLSD
[BER_JMLSD_SNR,~] = JMLSD_th( ...
    SNR_dB,a,M,T,channel,m_factor);

% SIC - Near user U1
BER_SIC_U1_SNR = SIC_nearUser_th( ...
    SNR_dB,a,M,T,channel,m_factor);

% SIC - Far user U2
BER_SIC_U2_SNR = SIC_farUser_th( ...
    SNR_dB,a,M,T,channel,m_factor);


%% ============================================================
% FIGURE 6(b): BER versus alpha1
% SNR = 15 dB
% ============================================================

SNR_alpha = 15;

alpha1 = 0.05:0.01:0.45;
alpha2 = 1-alpha1;

BER_JMLSD_ALPHA = zeros(2,length(alpha1));
BER_SIC_U1_ALPHA = zeros(1,length(alpha1));
BER_SIC_U2_ALPHA = zeros(1,length(alpha1));

for k = 1:length(alpha1)

    a_current = [alpha1(k) alpha2(k)];

    % JMLSD
    [temp_JMLSD,~] = JMLSD_th( ...
        SNR_alpha,a_current,M,T,channel,m_factor);

    BER_JMLSD_ALPHA(:,k) = temp_JMLSD(:);

    % SIC U1
    BER_SIC_U1_ALPHA(k) = SIC_nearUser_th( ...
        SNR_alpha,a_current,M,T,channel,m_factor);

    % SIC U2
    BER_SIC_U2_ALPHA(k) = SIC_farUser_th( ...
        SNR_alpha,a_current,M,T,channel,m_factor);

end


%% ============================================================
% SAVE
% ============================================================

if ~exist('results','dir')
    mkdir('results');
end
 
save(fullfile('results','results.mat'), ...
    'M','T','channel','m_factor', ...n
    'SNR_dB','a', ...
    'BER_JMLSD_SNR', ...
    'BER_SIC_U1_SNR', ...
    'BER_SIC_U2_SNR', ...
    'SNR_alpha','alpha1','alpha2', ...
    'BER_JMLSD_ALPHA', ...
    'BER_SIC_U1_ALPHA', ...
    'BER_SIC_U2_ALPHA');

disp('Figure 6 data generated successfully.');