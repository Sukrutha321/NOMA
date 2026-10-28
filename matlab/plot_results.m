clc;
clear;
close all;

load('results\results.mat');

%% ============================================================
% FIGURE 6
% ============================================================

figure('Color','white','Position',[100 100 900 450]);


%% ============================================================
% LEFT: BER vs SNR
% ============================================================

subplot(1,2,1);

hold on;

% Analytical curves
h1 = semilogy(SNR_dB,BER_SIC_U1_SNR, ...
    'b-','LineWidth',1.2);

h2 = semilogy(SNR_dB,BER_SIC_U2_SNR, ...
    'r-','LineWidth',1.2);

h3 = semilogy(SNR_dB,BER_JMLSD_SNR(1,:), ...
    'b--','LineWidth',1.2);

h4 = semilogy(SNR_dB,BER_JMLSD_SNR(2,:), ...
    'r:','LineWidth',1.2);


%% Markers
% These markers are placed on the analytical curves
% for visual comparison with the reference figure.

idx = [1 6 11 16];

semilogy(SNR_dB(idx),BER_SIC_U1_SNR(idx), ...
    'b*','MarkerSize',7,'LineStyle','none');

semilogy(SNR_dB(idx),BER_SIC_U2_SNR(idx), ...
    'r*','MarkerSize',7,'LineStyle','none');

semilogy(SNR_dB(idx),BER_JMLSD_SNR(1,idx), ...
    'bo','MarkerSize',6,'LineStyle','none');

semilogy(SNR_dB(idx),BER_JMLSD_SNR(2,idx), ...
    'ro','MarkerSize',6,'LineStyle','none');


%% Formatting

set(gca,'YScale','log');

xlim([0 20]);
ylim([1e-4 1]);

xlabel('SNR (dB)');
ylabel('P_{B_n}');

grid on;
grid minor;

set(gca,'FontSize',10);

title('(a)');


%% ============================================================
% RIGHT: BER vs ALPHA1
% ============================================================

subplot(1,2,2);

hold on;

% Analytical curves
semilogy(alpha1,BER_SIC_U1_ALPHA, ...
    'b-','LineWidth',1.2);

semilogy(alpha1,BER_SIC_U2_ALPHA, ...
    'r-','LineWidth',1.2);

semilogy(alpha1,BER_JMLSD_ALPHA(1,:), ...
    'b--','LineWidth',1.2);

semilogy(alpha1,BER_JMLSD_ALPHA(2,:), ...
    'r:','LineWidth',1.2);


%% Markers

idx2 = [1 6 11 16 21 26 31 36 41];

semilogy(alpha1(idx2),BER_SIC_U1_ALPHA(idx2), ...
    'b*','MarkerSize',7,'LineStyle','none');

semilogy(alpha1(idx2),BER_SIC_U2_ALPHA(idx2), ...
    'r*','MarkerSize',7,'LineStyle','none');

semilogy(alpha1(idx2),BER_JMLSD_ALPHA(1,idx2), ...
    'bo','MarkerSize',6,'LineStyle','none');

semilogy(alpha1(idx2),BER_JMLSD_ALPHA(2,idx2), ...
    'ro','MarkerSize',6,'LineStyle','none');


%% Formatting

set(gca,'YScale','log');

xlim([0 0.5]);
ylim([1e-4 1]);

xlabel('\alpha_1');
ylabel('P_{B_n}');

grid on;
grid minor;

set(gca,'FontSize',10);

title('(b)');


%% ============================================================
% LEGEND
% ============================================================

legend([h1 h2 h3 h4], ...
    'U_1: MR-SIC (42)', ...
    'U_2: MR-SIC (36)', ...
    'U_1: JMLSD (22)', ...
    'U_2: JMLSD (22)', ...
    'Location','northeast');