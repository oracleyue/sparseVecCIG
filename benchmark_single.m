% The script computes the BIC criterion to determine lambda, the
% regularization parameter.

% Copyright [2019] <oracleyue>
% Last modified on 24 Jun 2019


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

eTime = zeros(2, length(pList));
fprintf('Speed benchmark (elapsed CPU time):\n')
for k = 1:length(pList)

    % Data
    p = pList(k);
    dL = dLCell{k};
    Omega = sprandOm(dL, [.3 .8]);
    Sigma = inv(Omega);
    d = sum(dL); dList(k) = d;
    T = 10 * d;
    X = mvnrnd(zeros(T,d), Sigma);
    % sample covariance, normalized by T
    S = cov(X, 1);
    lambda = .5;

    fprintf('  [%2d]: #blocks=%2d, dim=%3d\n', k, p, d);
    aTimer = tic;
    [OmegaHat, ~] = bcdSpML(S, dL, lambda, [1e-3 10]);
    cgTime = toc(aTimer);
    fprintf('        CG   : %.6fs \n', cgTime);
    agTimer = tic;
    [F, Time, OmegaGoran] = Algorithm(speye(d), S, dL, lambda, 10, 1e-3, 0);
    goranTime = toc(agTimer);
    fprintf('        Goran: %.6fs \n', goranTime);

    eTime(:, k) = [cgTime; goranTime];
end
fprintf('End.\n')

%% visualization
fig_hl = figure;
subplot(1,2,1)
semilogy(pList, eTime(1, :), 'o--', pList, eTime(2, :), '^--');
legend('bcdpMLcg', 'Goran', 'location', 'northwest')
xlim([min(pList), max(pList)]);
xlabel('#blocks on diagonal');
ylabel('CPU time (s)')

subplot(1,2,2)
semilogy(dList, eTime(1, :), 'o--', dList, eTime(2, :), '^--');
legend('bcdpMLcg', 'Goran', 'location', 'northwest')
xlim([min(dList), max(dList)]);
xlabel('dimensions');
ylabel('CPU time (s)')

% save as pdf
set(fig_hl,'Units','Inches', 'Position', [3.5417 4.3021 8.7396 3.0625]);
pos = get(fig_hl,'Position');
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hl, 'benchmark_result', '-dpdf', '-r0')
