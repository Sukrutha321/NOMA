
%% Optimal Power Allocation for MR-NOMA
% Based on the optimization criterion in the reference paper.
% Requires:
%   SIC_nearUser_th.m
%   SIC_farUser_th.m
%   alpha1Max2UE.m
%
% IMPORTANT:
% Verify that SIC_nearUser_th supports M = [4 4], [16 4],
% and [64 4]. The available version was previously found to
% hardcode M = [2 2], so it must be corrected before using
% the QAM cases.

clc;
clear;
close all;

%% Paper parameters
zeta = 1e-2;
T2_values = [1 2 4 8];

M1_values = [4 16 64];
M2 = 4;
gamma1_values = [16 19 28];  % dB
gamma2_dB = 8;               % dB

channel = 'AWGN';
m_factor = [];

numCases = length(M1_values);
numT = length(T2_values);

alphaStar_U1 = nan(numCases,numT);
BER_U1_opt = nan(numCases,numT);
BER_U2_at_U1opt = nan(numCases,numT);

alphaStar_U2 = nan(numCases,numT);
BER_U2_opt = nan(numCases,numT);
BER_U1_at_U2opt = nan(numCases,numT);

%% Search resolution for alpha_1
Nalpha = 2000;

for c = 1:numCases

    M = [M1_values(c), M2];
    b = log2(M);

    % Maximum allocation allowed for reliable SIC
    alpha1_max = alpha1Max2UE(b);

    % Keep alpha_1 strictly inside the allowed interval
    alphaGrid = linspace(1e-4, ...
        min(alpha1_max - 1e-4, 0.999), Nalpha);

    for k = 1:numT

        T = [1 T2_values(k)];

        PB1 = nan(size(alphaGrid));
        PB2 = nan(size(alphaGrid));

        %% Evaluate both users' BER
        for j = 1:length(alphaGrid)

            a1 = alphaGrid(j);
            a = [a1, 1-a1];

            PB1(j) = SIC_nearUser_th( ...
                gamma1_values(c), a, M, T, channel, m_factor);

            PB2(j) = SIC_farUser_th( ...
                gamma2_dB, a, M, T, channel, m_factor);
        end

        %% Optimization 1:
        % Minimize near-user BER subject to PB2 <= zeta
        feasible1 = find(PB2 <= zeta);

        if ~isempty(feasible1)
            % PB1 decreases with alpha_1 in the paper's formulation
            idx = feasible1(end);

            alphaStar_U1(c,k) = alphaGrid(idx);
            BER_U1_opt(c,k) = PB1(idx);
            BER_U2_at_U1opt(c,k) = PB2(idx);
        end

        %% Optimization 2:
        % Minimize far-user BER subject to PB1 <= zeta
        feasible2 = find(PB1 <= zeta);

        if ~isempty(feasible2)
            % PB2 increases with alpha_1 in the paper's formulation
            idx = feasible2(1);

            alphaStar_U2(c,k) = alphaGrid(idx);
            BER_U2_opt(c,k) = PB2(idx);
            BER_U1_at_U2opt(c,k) = PB1(idx);
        end

    end
end

%% Plot 1: Optimal alpha_1 versus T_2
figure('Color','w');
hold on;

for c = 1:numCases
    plot(T2_values, alphaStar_U1(c,:), '-o', ...
        'LineWidth', 1.7, ...
        'DisplayName', sprintf('M = [%d,4]', M1_values(c)));
end

xlabel('Far-user symbol duration T_2');
ylabel('Optimal power allocation \alpha_1^*');
title('Optimal \alpha_1^* versus T_2');
legend('Location','best');
grid on;
xticks(T2_values);

%% Plot 2: BER at the near-user-optimal allocation
figure('Color','w');
hold on;

for c = 1:numCases
    semilogy(T2_values, max(BER_U1_opt(c,:), eps), '-o', ...
        'LineWidth', 1.7, ...
        'DisplayName', sprintf('Near user, M = [%d,4]', M1_values(c)));

    semilogy(T2_values, max(BER_U2_at_U1opt(c,:), eps), '--s', ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Far user constraint, M = [%d,4]', M1_values(c)));
end

yline(zeta, 'k:', 'BER threshold \zeta');
xlabel('Far-user symbol duration T_2');
ylabel('BER');
title('BER at the near-user-optimal allocation');
legend('Location','best');
grid on;
xticks(T2_values);

%% Save results
if ~exist('results','dir')
    mkdir('results');
end

save(fullfile('results','optimal_time_assignment.mat'), ...
    'T2_values','M1_values','gamma1_values','gamma2_dB','zeta', ...
    'alphaStar_U1','BER_U1_opt','BER_U2_at_U1opt', ...
    'alphaStar_U2','BER_U2_opt','BER_U1_at_U2opt');

disp('Near-user optimization:');
disp(array2table(alphaStar_U1, ...
    'VariableNames', {'T2_1','T2_2','T2_4','T2_8'}, ...
    'RowNames', {'M1_4','M1_16','M1_64'}));

disp('Results saved in the results folder.');
