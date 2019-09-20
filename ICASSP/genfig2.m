% Require to run "benchmark.m" first to get the .mat result file.

% Copyright [2019] <oracleyue>
% Last modified on 14 Aug 2019


close all;

% load data
load('./bm_results.mat');

% choose x-axis
xType = 'dim';
% xType = 'blk';

% plotting
fig_hl = figure;
cmap = colormap('lines');
grayZ = cmap(1,:);
grayG = cmap(2,:);
alpha = .08;

datTimeZ = eTimeMatrix(1:numL, :);
datTimeG = eTimeMatrix(numL+1:2*numL, :);
datTimeZ = log10(datTimeZ);
datTimeG = log10(datTimeG);
eTime = log10(eTime);
meanTimeZ = mean(datTimeZ);
stdTimeZ  = std(datTimeZ);
meanTimeG = mean(datTimeG);
stdTimeG  = std(datTimeG);

switch xType
  case 'blk'    % x-axis: #diagonal blocks
    plot(pList, (meanTimeZ), 'o--', 'Color', cmap(1,:), 'LineWidth', 1.2);
    hold on
    plot(pList, (meanTimeG), 's--', 'Color', cmap(2,:), 'LineWidth', 1.2);
    plot(pList, (eTime(1,:)), '*', 'Color', cmap(1,:));
    plot(pList, (eTime(2,:)), '*', 'Color', cmap(2,:));
    xlim([min(pList), max(pList)]);
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
    plot(dList, (meanTimeZ), 'o--', 'Color', cmap(1,:), 'LineWidth', 1.2);
    hold on
    plot(dList, (meanTimeG), 's--', 'Color', cmap(2,:), 'LineWidth', 1.2);
    plot(dList, (eTime(1,:)), '*', 'Color', cmap(1,:));
    plot(dList, (eTime(2,:)), '*', 'Color', cmap(2,:));
    xlim([min(dList), max(dList)]);
    xlabel('dimensions');
    ylabel('CPU time in log10 (s)');

    % draw deviation regionlog10
    polyX = [dList fliplr(dList)];
    polyYZ = [meanTimeZ-stdTimeZ fliplr(meanTimeZ+stdTimeZ)];
    ylim([-3 3]); polyYZ(1) = -3; % since polyYZ is log10 of a negative value
    patchZ = fill(polyX, polyYZ, grayZ);
    set(patchZ, 'edgecolor', 'none');
    set(patchZ, 'FaceAlpha', alpha);
    polyYG = [meanTimeG-stdTimeG fliplr(meanTimeG+stdTimeG)];
    patchG = fill(polyX, polyYG, grayG);
    set(patchG, 'edgecolor', 'none');
    set(patchG, 'FaceAlpha', alpha);

    leg_hl = legend('MLCG', 'Goran', 'MLCG (best)', ...
                    'Goran (best)', 'location', 'northwest');
    legstr = get(leg_hl, 'string');
    set(leg_hl, 'string', legstr(1:4))
end

%% save as pdf
pos = [3.5417 4.3021 4 2.5];
set(fig_hl,'Units','Inches', 'Position', pos);
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hl, 'bmk-plot', '-dpdf', '-r0')
