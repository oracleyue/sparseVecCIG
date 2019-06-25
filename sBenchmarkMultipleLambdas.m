% The script is to benchmark speeds of two implementations (CG-embeded
% vs. Goran's), by running the algorithms for a range of lambda. It
% provides the statistical plot of time costs over a range of lambda and
% the plot of time cost with the BIC-optimal lambda.
%
% Notes: It may not be reasonable to benchmark time costs by running with an
% arbitrary fixed lambda, since the choice of lambda might effect the
% convergence of algorithms and hence time costs.

% Copyright [2019] <oracleyue>
% Last modified on 25 Jun 2019


clear all; close all;

% Search paths
addpath('./extern'); % external libraries
addpath('./Goran');  % using Goran's
addpath('./tools');

% Init seed
rng(2);

%% Benchmark for different dimensions
pList = 4:2:20;
dList = zeros(size(pList));
dLCell = cell(size(pList));
for k = 1:length(pList)
    if k == 1
        dLCell{k} = randi(5, pList(1), 1)*3;
    else
        dLCell{k} = [dLCell{k-1}; randi(5, pList(k)-pList(k-1), 1)*3];
    end
end

% set lambda range
numL = 40;
lambdaList = logspace(-2, 0, numL);

% init saving variables
lambdaBest = zeros(2, length(pList));     % 1: 'zyue'; 2: 'goran'
infoCrit = zeros(2*numL, length(pList));  % 1:numL: 'zyue';
                                          % numL+1:2*numL: 'goran'
eTimeMatrix = zeros(2*numL, length(pList));
eTime = zeros(2, length(pList));

fprintf('Speed benchmark (elapsed CPU time):\n')
for k = 1:length(pList)
    % data
    p = pList(k);
    dL = dLCell{k};
    Omega = sprandOm(dL, [.3 .8]);
    Sigma = inv(Omega);
    d = sum(dL); dList(k) = d;
    N = 10 * d;
    X = mvnrnd(zeros(N,d), Sigma);
    % sample covariance, normalized by T
    S = cov(X, 1);

    fprintf('  [%2d]: #blocks=%2d, dim=%3d\n', k, p, d);

    % ML with different lambdas
    [lambdaZ, ICsZ, eTimeZ] = calcLambda(S, dL, N, lambdaList, 'AIC', 'zyue');
    fprintf('        CG   : %.6fs \n', sum(eTimeZ));
    [lambdaG, ICsG, eTimeG] = calcLambda(S, dL, N, lambdaList, 'AIC', 'goran');
    fprintf('        Goran: %.6fs \n', sum(eTimeG));

    lambdaBest(:,k) = [lambdaZ; lambdaG];
    infoCrit(1:numL,k) = ICsZ;
    infoCrit(numL+1:2*numL,k) = ICsG;
    eTimeMatrix(1:numL, k) = eTimeZ;
    eTimeMatrix(numL+1:2*numL, k) = eTimeG;
    eTime(1, k) = eTimeZ(find(lambdaList == lambdaZ));
    eTime(2, k) = eTimeG(find(lambdaList == lambdaG));
end
fprintf('End.\n')
save('benchmarkMulti_results.mat');

%% visualization
fig_hl = figure;
cmap = colormap('lines');
grayZ = cmap(1,:);
grayG = cmap(2,:);
% gmap = colormap('gray');
% grayZ = gmap(30, :); grayG = grayZ;
alpha = .08;

meanTimeZ = mean(log10(eTimeMatrix(1:numL, :)));
stdTimeZ = std(log10(eTimeMatrix(1:numL, :)));
meanTimeG = mean(log10(eTimeMatrix(numL+1:2*numL, :)));
stdTimeG = std(log10(eTimeMatrix(numL+1:2*numL, :)));

subplot(1,2,1)
plot(pList, meanTimeZ, 'o--', 'Color', cmap(1,:), 'LineWidth', 1.2);
hold on
plot(pList, meanTimeG, 's--', 'Color', cmap(2,:), 'LineWidth', 1.2);
plot(pList, log10(eTime(1,:)), '*', 'Color', cmap(1,:));
plot(pList, log10(eTime(2,:)), '*', 'Color', cmap(2,:));
xlim([min(pList), max(pList)]);
xlabel('#blocks on diagonal');
ylabel('Time in log scale (s)')

polyX = [pList fliplr(pList)];
polyYZ = [meanTimeZ-stdTimeZ fliplr(meanTimeZ+stdTimeZ)];
patchZ = fill(polyX, polyYZ, grayZ);
set(patchZ, 'edgecolor', 'none');
set(patchZ, 'FaceAlpha', alpha);
polyYG = [meanTimeG-stdTimeG fliplr(meanTimeG+stdTimeG)];
patchG = fill(polyX, polyYG, grayG);
set(patchG, 'edgecolor', 'none');
set(patchG, 'FaceAlpha', alpha);

leg_hl = legend('MLCG', 'Goran', 'MLCG (best)', ...
                'Goran (best)', 'location', 'northwest')
legstr = get(leg_hl, 'string');
set(leg_hl, 'string', legstr(1:4))

subplot(1,2,2)
plot(dList, meanTimeZ, 'o--', 'Color', cmap(1,:), 'LineWidth', 1.2);
hold on
plot(dList, meanTimeG, 's--', 'Color', cmap(2,:), 'LineWidth', 1.2);
plot(dList, log10(eTime(1,:)), '*', 'Color', cmap(1,:));
plot(dList, log10(eTime(2,:)), '*', 'Color', cmap(2,:));
xlim([min(dList), max(dList)]);
xlabel('dimensions');
ylabel('Time in log scale (s)')

polyX = [dList fliplr(dList)];
polyYZ = [meanTimeZ-stdTimeZ fliplr(meanTimeZ+stdTimeZ)];
patchZ = fill(polyX, polyYZ, grayZ);
set(patchZ, 'edgecolor', 'none');
set(patchZ, 'FaceAlpha', alpha);
polyYG = [meanTimeG-stdTimeG fliplr(meanTimeG+stdTimeG)];
patchG = fill(polyX, polyYG, grayG);
set(patchG, 'edgecolor', 'none');
set(patchG, 'FaceAlpha', alpha);

leg_hl = legend('MLCG', 'Goran', 'MLCG (best)', ...
                'Goran (best)', 'location', 'northwest')
legstr = get(leg_hl, 'string');
set(leg_hl, 'string', legstr(1:4))

% save as pdf
set(fig_hl,'Units','Inches', 'Position', [3.5417 4.3021 8.7396 3.0625]);
pos = get(fig_hl,'Position');
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hl, 'benchmarkMulti_results', '-dpdf', '-r0')
