% Plot demo results

% Copyright [2019] <oracleyue>
% Last modified on 25 Sep 2019


clear all; close all;
load('./Goran_Demo.mat');

% visualization
fig_hl = figure;
set(gcf,'color','white');

[ha, pos] = tight_subplot(1,4, [.001 .005], [.01 .01], [.01 .01])

axes(ha(1));
imshowOm(Omega, 'raw');
title('Ground truth $\Omega$', ...
      'FontSize', 12, 'Interpreter','latex');

axes(ha(2));
imshowOm(OmHatL0, 'raw');
title('Estimation $\mathbf{\widehat{\Omega}}$', ...
      'FontSize', 12, 'Interpreter','latex');

axes(ha(3));
imshowOm(OmHatL0g, 'raw');
title('Estimation $\mathbf{\widehat{\Omega}}_{l_0}$', ...
      'FontSize', 12, 'Interpreter','latex');

axes(ha(4));
imshowOm(OmHatL1g, 'raw');
title('Estimation $\mathbf{\widehat{\Omega}}_{l_1}$', ...
      'FontSize', 12, 'Interpreter','latex');

% export as pdf
filename = ['./goran-demo', '.pdf'];
pos = [4.8542 5.8958 6.8958 2.2396];
set(fig_hl,'Units','Inches', 'Position', pos);
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hl, filename, '-dpdf', '-r0')
