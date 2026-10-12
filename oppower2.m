% oppower2.m
% Fig. 11-style optimal power allocation plots for two-user MR-NOMA.
% Requires SIC_nearUser_th.m, SIC_farUser_th.m, alpha1Max2UE.m and their helpers.
%
% Paper settings: zeta=1e-2, M2=4, M1=[4 16 64],
% gamma1=[16 19 28] dB, gamma2=8 dB.
%
% IMPORTANT: This calls the near-user BER routine at gamma1 and the far-user
% BER routine at gamma2. Verify that these routines interpret their first
% argument as the SNR of the user being evaluated. If both functions assume
% a common system SNR, this is not yet an exact reproduction of Fig. 11.
%
% Eq. (46): minimize PB1 subject to PB2 <= zeta and alpha1 < alpha1,max.
% Eq. (47): minimize PB2 subject to PB1 <= zeta and alpha1 < alpha1,max.
% Marker: square [4,4], circle [16,4], star [64,4].
clear; clc; close all;

M1_values = [4 16 64];
M2 = 4;
gamma1_dB = [16 19 28];
gamma2_dB = 8;
zeta = 1e-2;
T2_values = 1:12;          % Change to match the paper's displayed T2 values.
nGrid = 4001;
markers = {'s','o','*'};
pairLabels = {'[4,4]','[16,4]','[64,4]'};
nM = numel(M1_values);
nT = numel(T2_values);

alphaStar = nan(nM,nT,2); % third dimension: 1=Eq.(46), 2=Eq.(47)
PB1_at_star = nan(nM,nT,2);
PB2_at_star = nan(nM,nT,2);
alphaMaxAll = nan(nM,1);

for im = 1:nM
    M = [M1_values(im), M2];
    b = log2(M);
    alphaMax = alpha1Max2UE(b);
    if ~isscalar(alphaMax) || ~isfinite(alphaMax) || alphaMax <= 0
        error('Invalid alpha1,max for M=[%d,%d].',M(1),M(2));
    end
    alphaMaxAll(im) = alphaMax;
    epsA = max(1e-5,1e-4*alphaMax);
    aLo = epsA;
    aHi = alphaMax-epsA;
    if aHi <= aLo, error('Invalid alpha search interval.'); end
    alphaGrid = linspace(aLo,aHi,nGrid);

    for it = 1:nT
        T = [1 T2_values(it)];
        PB1 = nan(size(alphaGrid));
        PB2 = nan(size(alphaGrid));
        for ia = 1:numel(alphaGrid)
            a1 = alphaGrid(ia);
            a = [a1,1-a1];
            PB1(ia) = SIC_nearUser_th(gamma1_dB(im),a,M,T,'AWGN',[]);
            PB2(ia) = SIC_farUser_th(gamma2_dB,a,M,T,'AWGN',[]);
        end
        PB1(~isfinite(PB1) | PB1<0) = NaN;
        PB2(~isfinite(PB2) | PB2<0) = NaN;

        % Eq. (46): minimize PB1 with PB2 <= zeta.
        idxFeasible = find(PB2<=zeta & isfinite(PB1));
        if ~isempty(idxFeasible)
            [~,j] = min(PB1(idxFeasible));
            idx = idxFeasible(j);
            alphaStar(im,it,1) = alphaGrid(idx);
            PB1_at_star(im,it,1) = PB1(idx);
            PB2_at_star(im,it,1) = PB2(idx);
        end

        % Eq. (47): minimize PB2 with PB1 <= zeta.
        idxFeasible = find(PB1<=zeta & isfinite(PB2));
        if ~isempty(idxFeasible)
            [~,j] = min(PB2(idxFeasible));
            idx = idxFeasible(j);
            alphaStar(im,it,2) = alphaGrid(idx);
            PB1_at_star(im,it,2) = PB1(idx);
            PB2_at_star(im,it,2) = PB2(idx);
        end
    end
    fprintf('Finished [%d,%d]; alpha1,max = %.6g\\n',M(1),M(2),alphaMax);
end

% Figure 1: T2 versus alpha1* for both optimizations.
figure('Color','w','Name','Fig. 11(a) - Optimal alpha1');
hold on; box on; grid on;
for im = 1:nM
    plot(T2_values,squeeze(alphaStar(im,:,1)),['-' markers{im}], ...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName',[pairLabels{im} ' Eq. (46)']);
    plot(T2_values,squeeze(alphaStar(im,:,2)),['--' markers{im}], ...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName',[pairLabels{im} ' Eq. (47)']);
end
xlabel('T_2'); ylabel('\alpha_1^*');
title('Optimal power allocation versus T_2');
legend('Location','bestoutside'); set(gca,'FontSize',11);

% Figure 2: BER at alpha1*. Separate panels keep all curves readable.
figure('Color','w','Name','Fig. 11(b) - BER at optimal alpha');
subplot(2,1,1); hold on; box on; grid on;
for im = 1:nM
    semilogy(T2_values,squeeze(PB1_at_star(im,:,1)),['-' markers{im}], ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName',[pairLabels{im} ' P_{B1}, Eq. (46)']);
    semilogy(T2_values,squeeze(PB1_at_star(im,:,2)),[':' markers{im}], ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName',[pairLabels{im} ' P_{B1}, Eq. (47)']);
end
yline(zeta,'k--','\zeta = 10^{-2}');
xlabel('T_2'); ylabel('P_{B1}(\alpha_1^*)');
title('Near-user BER at the optimum'); legend('Location','eastoutside');
set(gca,'YScale','log','FontSize',9);

subplot(2,1,2); hold on; box on; grid on;
for im = 1:nM
    semilogy(T2_values,squeeze(PB2_at_star(im,:,1)),['-' markers{im}], ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName',[pairLabels{im} ' P_{B2}, Eq. (46)']);
    semilogy(T2_values,squeeze(PB2_at_star(im,:,2)),[':' markers{im}], ...
        'LineWidth',1.4,'MarkerSize',6,'DisplayName',[pairLabels{im} ' P_{B2}, Eq. (47)']);
end
yline(zeta,'k--','\zeta = 10^{-2}');
xlabel('T_2'); ylabel('P_{B2}(\alpha_1^*)');
title('Far-user BER at the optimum'); legend('Location','eastoutside');
set(gca,'YScale','log','FontSize',9);

results = struct('M1_values',M1_values,'M2',M2,'gamma1_dB',gamma1_dB, ...
    'gamma2_dB',gamma2_dB,'zeta',zeta,'T2_values',T2_values, ...
    'alphaMaxAll',alphaMaxAll,'alphaStar',alphaStar, ...
    'PB1_at_star',PB1_at_star,'PB2_at_star',PB2_at_star);
save('oppower2_results.mat','results');
disp('Finished. Numerical results saved to oppower2_results.mat.');
disp('NaN means no feasible alpha1 was found on the selected grid.');
