function plot_res_fig9_10(M_pairs)

clc;
close all;

%% =========================================================
% FIGURE 9
% Far-user MR-SIC BER versus SNR
%
% Example:
% plot_res_fig9_10([2 2; 4 4])
%
% or
% plot_res_fig9_10([16 16; 64 64])
% ==========================================================

%% SNR range
SNR_dB = -5:1:35;

%% Channel
channel = 'AWGN';

%% Additional parameter
m_factor = [];

%% Fixed symbol-duration cases
T_cases = {
    [1 1]
    [1 2]
    [1 4]
    [1 8]
};

%% ---------------------------------------------------------
% Check input
%% ---------------------------------------------------------

if nargin < 1
    error('Example: plot_res_fig9_10([2 2; 4 4])');
end

if size(M_pairs,1) ~= 2 || size(M_pairs,2) ~= 2
    error('M_pairs must be a 2x2 matrix.');
end

%% ---------------------------------------------------------
% Plot styles
%
% Line style -> M pair
% Marker     -> T case
%% ---------------------------------------------------------

lineStyles = {'-','--'};

markers = {'s','o','*','^'};
% Colors correspond to T cases
T_colors = {'b','r','c','m'};
%% ---------------------------------------------------------
% Create ONE figure
%% ---------------------------------------------------------

figure( ...
    'Color','white', ...
    'Position',[100 100 900 650]);

hold on;

%% ---------------------------------------------------------
% Calculate and plot all 8 curves
%% ---------------------------------------------------------

for m = 1:2

    %% Current modulation pair
    M = M_pairs(m,:);

    %% IMPORTANT:
    % Power allocation is calculated from M
    % exactly as in the paper's code
    a = equal_distance(log2(M)).';

    for t = 1:4

        %% Current symbol durations
        T = T_cases{t};

        %% Far-user MR-SIC analytical BER
        BER_U2 = SIC_farUser_th( ...
            SNR_dB, ...
            a, ...
            M, ...
            T, ...
            channel, ...
            m_factor);

        %% Protect logarithmic plot
        BER_U2(BER_U2 <= 0) = 1e-10;

        %% Plot
        semilogy( ...
    SNR_dB, ...
    BER_U2, ...
    'Color',T_colors{t}, ...
    'LineStyle',lineStyles{m}, ...
    'Marker',markers{t}, ...
    'LineWidth',1.2, ...
    'MarkerSize',5);

    end

end

%% ---------------------------------------------------------
% Axis formatting
%% ---------------------------------------------------------

set(gca,'YScale','log');

xlim([-5 20]);
ylim([1e-4 1]);

xlabel('SNR (dB)');
ylabel('P_{B_2}');

grid on;
grid minor;

set(gca,'FontSize',11);

%% ---------------------------------------------------------
% Legend
%% ---------------------------------------------------------



%% ---------------------------------------------------------
% Legend — paper-style colors
%% ---------------------------------------------------------

hLine1 = plot(nan,nan,'k-', ...
    'LineWidth',1.2);

hLine2 = plot(nan,nan,'k--', ...
    'LineWidth',1.2);

hT1 = plot(nan,nan,'bs', ...
    'LineStyle','none','MarkerSize',5);

hT2 = plot(nan,nan,'ro', ...
    'LineStyle','none','MarkerSize',6);

hT4 = plot(nan,nan,'c*', ...
    'LineStyle','none','MarkerSize',6);

hT8 = plot(nan,nan,'m^', ...
    'LineStyle','none','MarkerSize',6);

legend( ...
    [hLine1 hLine2 hT1 hT2 hT4 hT8], ...
    { ...
    sprintf('[M_1,M_2] = [%d,%d]',M_pairs(1,1),M_pairs(1,2)), ...
    sprintf('[M_1,M_2] = [%d,%d]',M_pairs(2,1),M_pairs(2,2)), ...
    '[T_1,T_2] = [1,1]', ...
    '[T_1,T_2] = [1,2]', ...
    '[T_1,T_2] = [1,4]', ...
    '[T_1,T_2] = [1,8]' ...
    }, ...
    'Location','southwest');
%% ---------------------------------------------------------
% Title
%% ---------------------------------------------------------

title('Far-User MR-SIC BER vs SNR');

end