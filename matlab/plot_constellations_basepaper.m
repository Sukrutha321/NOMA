clc;
clear;
close all;

M_pairs = [4 4; 16 4; 64 4];

% Muted paper-style red
pointColor = [0.78 0.24 0.22];

figure('Color','w','Position',[80 80 1000 1050]);

% Generate U2 (4-QAM) once
[s2,~,~] = rectGrayQAM_mod(4);
s2 = s2 ./ sqrt(mean(abs(s2).^2));

% Use equal_distance for [M1,4] only once per row
subplot(3,2,2); % fixed column for U2
scatter(real(s2),imag(s2),32,'filled', ...
    'MarkerFaceColor',pointColor,'MarkerEdgeColor',pointColor);
axis equal; grid on; grid minor;
xlabel('In-phase'); ylabel('Quadrature');
title('U_2: 4-QAM (scaled)');
set(gca,'FontSize',9);

for row = 1:size(M_pairs,1)
    M = M_pairs(row,:);

    % Allocation coefficients
    a = equal_distance(log2(M));
    a = a(:).';
    a = a ./ sum(a);

    % Generate U1 constellation
    [s1,~,~] = rectGrayQAM_mod(M(1));
    s1 = s1 ./ sqrt(mean(abs(s1).^2));

    % Apply scaling
    s1_scaled = sqrt(a(1)) .* s1;
    s2_scaled = sqrt(a(2)) .* s2;

    % Composite NOMA constellation
    [S1,S2] = ndgrid(s1,s2);
    sNOMA = sqrt(a(1)).*S1(:) + sqrt(a(2)).*S2(:);

    % Plot U1
    subplot(3,2,(row-1)*2+1);
    scatter(real(s1_scaled),imag(s1_scaled),32,'filled', ...
        'MarkerFaceColor',pointColor,'MarkerEdgeColor',pointColor);
    axis equal; grid on; grid minor;
    xlabel('In-phase'); ylabel('Quadrature');
    title(sprintf('U_1: %d-QAM (scaled)',M(1)));
    set(gca,'FontSize',9);

    % Plot NOMA superposition
    subplot(3,2,(row-1)*2+2);
    scatter(real(sNOMA),imag(sNOMA),18,'filled', ...
        'MarkerFaceColor',pointColor,'MarkerEdgeColor',pointColor);
    axis equal; grid on; grid minor;
    xlabel('In-phase'); ylabel('Quadrature');
    title(sprintf('NOMA superposition [%d,%d]',M(1),M(2)));
    set(gca,'FontSize',9);

    fprintf('M = [%d %d], alpha1 = %.6f, alpha2 = %.6f\n', ...
        M(1),M(2),a(1),a(2));
end

sgtitle('MR-NOMA Constellations (U_2 shown once)');
