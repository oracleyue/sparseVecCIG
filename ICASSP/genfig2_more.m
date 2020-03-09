% Require to run "benchmark.m" first to get the .mat result file.

% Copyright [2019] <oracleyue>
% Last modified on 14 Aug 2019


close all; clear all;

% load data
load('./bm_results.mat');
load('./bm_results_zyue.mat');

% choose x-axis
xType = 'dim';
% xType = 'blk';

% plotting
fig_hl = figure;
cmap = colormap('lines');
grayZ = cmap(1,:);
grayG = cmap(2,:);
alpha = .08;

datTimeZ = eTimeMatrixZ;
datTimeG = eTimeMatrix(numL+1:2*numL, :);
datTimeZ = log10(datTimeZ);
datTimeG = log10(datTimeG);
eTimeZ = log10(eTimeZ);
eTimeG = log10(eTime(2,:));
meanTimeZ = mean(datTimeZ);
stdTimeZ  = std(datTimeZ);
meanTimeG = mean(datTimeG);
stdTimeG  = std(datTimeG);

switch xType
  case 'blk'    % x-axis: #diagonal blocks
    plot(pListZ, meanTimeZ, 'o--', 'Color', cmap(1,:), 'LineWidth', 1.2);
    hold on
    plot(pList,  meanTimeG, 's--', 'Color', cmap(2,:), 'LineWidth', 1.2);
    plot(pListZ, eTimeZ, '*', 'Color', cmap(1,:));
    plot(pList,  eTimeG, '*', 'Color', cmap(2,:));
    xlim([min(pListZ), max(pListZ)]);
    xlabel('#blocks on diagonal');
    ylabel('CPU time in log10 (s)');

    % draw deviation region
    polyX = [pList fliplr(pList)];
    polyYZ = [meanTimeZ-stdTimeZ fliplr(meanTimeZ+stdTimeZ)];
    patchZ = fill(polyX, polyYZ, grayZ);
    set(patchZ, 'edgecolor', 'none');
    set(patchZ, 'FaceAlpha', alpha);
    polyYG = [meanTimeG-stdTimeG fliplr(meanTimeG+stdTimeG)];
    patchG = fill(polyX, polyYG, grayG);
    set(patchG, 'edgecolor', 'none');log10
    set(patchG, 'FaceAlpha', alpha);

    leg_hl = legend('MLCG', 'Goran', 'MLCG (best)', ...
                    'Goran (best)', 'location', 'northwest');
    legstr = get(leg_hl, 'string');
    set(leg_hl, 'string', legstr(1:4))

  case 'dim'    % x-axis: dimension of X
    plot(dListZ, meanTimeZ, 'o--', 'Color', cmap(1,:), 'LineWidth', 1.2);
    hold on
    plot(dList,  meanTimeG, 's--', 'Color', cmap(2,:), 'LineWidth', 1.2);
    plot(dListZ, eTimeZ, '*', 'Color', cmap(1,:));
    plot(dList,  eTimeG, '*', 'Color', cmap(2,:));
    xlim([0, max(dListZ)]);
    xticks([0:50:200 300:100:800 1200 1700 2200]);
    xticklabels({'0', '', '', '', '200', '', '400', '', '600', '', '800', '1200', '1700', '2200'});
    xlabel('total dimension $d$', 'interpreter','latex');
    ylabel('CPU time in $\log_{10}$ (s)', 'interpreter','latex');
    grid on
    
    % draw deviation regionlog10
    polyXZ = [dListZ fliplr(dListZ)];
    polyYZ = [meanTimeZ-stdTimeZ fliplr(meanTimeZ+stdTimeZ)];
    patchZ = fill(polyXZ, polyYZ, grayZ);
    set(patchZ, 'edgecolor', 'none');
    set(patchZ, 'FaceAlpha', alpha);
    polyXG = [dList fliplr(dList)];
    polyYG = [meanTimeG-stdTimeG fliplr(meanTimeG+stdTimeG)];
    patchG = fill(polyXG, polyYG, grayG);
    set(patchG, 'edgecolor', 'none');
    set(patchG, 'FaceAlpha', alpha);

    leg_hl = legend('FB$l_0$-SIC', 'B$l_0$-SIC', 'FB$l_0$-SIC (best)', ...
                    'B$l_0$-SIC (best)', 'location', 'southeast', 'interpreter','latex');
    legstr = get(leg_hl, 'string');
    set(leg_hl, 'string', legstr(1:4))
end

%% save as pdf
pos = [3.5417 4.3021 4 2.5];
set(fig_hl,'Units','Inches', 'Position', pos);
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hl, 'bmk-plot', '-dpdf', '-r0')
