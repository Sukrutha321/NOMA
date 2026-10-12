%% Optimal power allocation for MR-NOMA  (problems (43) and (44) of the paper)
% Fixes the earlier draft: the near-user BER used for the optimization is the
% perfect-SIC lower bound  P_B1^lb(alpha1) = singleUser_th(gamma1, alpha1, M1)
% (Eq. (45)), NOT SIC_nearUser_th (which is hard-coded to M = [2 2]).
% Requires: SIC_farUser_th.m, singleUser_th.m, alpha1Max2UE.m,
%           nchoosek_sVer2.m, rectGrayQAM_mod.m
clc; clear; close all;

zeta       = 1e-2;
T2_values  = [1 2 4 8];
M1_values  = [4 16 64];  M2 = 4;
gamma1_dB  = [16 19 28]; gamma2_dB = 8;
channel    = 'AWGN';     m_factor = [];
Nalpha     = 2000;

nC = numel(M1_values); nT = numel(T2_values);
a43 = nan(nC,nT); PB1_43 = nan(nC,nT); PB2_43 = nan(nC,nT);   % problem (43)
a44 = nan(nC,nT); PB1_44 = nan(nC,nT); PB2_44 = nan(nC,nT);   % problem (44)

for c = 1:nC
    M = [M1_values(c) M2];
    amax = alpha1Max2UE(log2(M));
    g = linspace(1e-4, min(amax-1e-4,0.999), Nalpha);

    % near-user BER: perfect-SIC lower bound, independent of T2
    PB1 = arrayfun(@(x) singleUser_th(gamma1_dB(c), x, M(1), channel, m_factor), g);

    for k = 1:nT
        T = [1 T2_values(k)];
        PB2 = arrayfun(@(x) SIC_farUser_th(gamma2_dB,[x 1-x],M,T,channel,m_factor), g);

        % (43): minimize PB1 s.t. PB2 <= zeta  -> largest feasible alpha1
        i1 = find(PB2 <= zeta, 1, 'last');
        if ~isempty(i1)
            a43(c,k)=g(i1); PB1_43(c,k)=PB1(i1); PB2_43(c,k)=PB2(i1);
        end
        % (44): minimize PB2 s.t. PB1 <= zeta  -> smallest feasible alpha1
        i2 = find(PB1 <= zeta, 1, 'first');
        if ~isempty(i2)
            a44(c,k)=g(i2); PB1_44(c,k)=PB1(i2); PB2_44(c,k)=PB2(i2);
        end
    end
end

%% Plot style: solid = (43), dashed = (44); square [4,4], circle [16,4], star [64,4]
mk  = {'s','o','*'};  col = {[0 .447 .741],[.851 .325 .098],[.18 .55 .18]};
lbl = {'[4,4]','[16,4]','[64,4]'};

%% Figure A: alpha1* versus T2
figure('Color','w','Position',[100 100 560 430]); hold on;
for c = 1:nC
    plot(T2_values,a43(c,:),'-' ,'Marker',mk{c},'Color',col{c},'LineWidth',1.3,'MarkerSize',7);
    plot(T2_values,a44(c,:),'--','Marker',mk{c},'Color',col{c},'LineWidth',1.3,'MarkerSize',7);
end
set(gca,'XScale','log','XTick',T2_values,'XTickLabel',T2_values,'FontSize',11);
xlabel('T_2'); ylabel('\alpha_1^*'); grid on; ylim([0 0.55]);
h = [plot(nan,nan,'k-'), plot(nan,nan,'k--')];
for c = 1:nC, h(end+1) = plot(nan,nan,'LineStyle','none','Marker',mk{c},'Color',col{c},'MarkerSize',7); end
legend(h,[{'Problem (43)','Problem (44)'},lbl],'Location','southeast');

%% Figure B: P_B1(alpha1*) and P_B2(alpha1*) versus T2
figure('Color','w','Position',[700 100 560 470]); hold on;
cB1 = [0 .447 .741]; cB2 = [.753 .224 .169];
for c = 1:nC
    semilogy(T2_values,max(PB1_43(c,:),eps),'-' ,'Marker',mk{c},'Color',cB1,'LineWidth',1.3,'MarkerSize',7);
    semilogy(T2_values,max(PB2_43(c,:),eps),'-' ,'Marker',mk{c},'Color',cB2,'LineWidth',1.3,'MarkerSize',7);
    semilogy(T2_values,max(PB1_44(c,:),eps),'--','Marker',mk{c},'Color',cB1,'LineWidth',1.3,'MarkerSize',7);
    semilogy(T2_values,max(PB2_44(c,:),eps),'--','Marker',mk{c},'Color',cB2,'LineWidth',1.3,'MarkerSize',7);
end
yline(zeta,'k:','\zeta = 10^{-2}');
set(gca,'YScale','log','XScale','log','XTick',T2_values,'XTickLabel',T2_values,'FontSize',11);
xlabel('T_2'); ylabel('P_{B_n}(\alpha_1^*)'); grid on; ylim([1e-7 1]);
h = [plot(nan,nan,'-','Color',cB1), plot(nan,nan,'-','Color',cB2), plot(nan,nan,'k-'), plot(nan,nan,'k--')];
for c = 1:nC, h(end+1) = plot(nan,nan,'LineStyle','none','Marker',mk{c},'Color',[.4 .4 .4],'MarkerSize',7); end
legend(h,[{'P_{B_1}(\alpha_1^*)','P_{B_2}(\alpha_1^*)','Problem (43)','Problem (44)'},lbl],'Location','southwest','NumColumns',2);

if ~exist('results','dir'), mkdir('results'); end
save(fullfile('results','optimal_final.mat'),'T2_values','M1_values','zeta','a43','PB1_43','PB2_43','a44','PB1_44','PB2_44');
