n_homo = [5:15];

t_homo = [
 5.17 
 12.37 
 25.80 
 50.48 
 90.78 
 161.54 
 263.77 
 411.40 
 634.18 
 1136.59 
 1451.06 
];

n_mep = [5:15]';

t_mep = [
 0.10 
 0.43 
 1.55 
 5.97 
 18.73 
 47.00 
 112.62 
 235.92 
 473.54 
 986.88 
 1762.43 
 ];

% same size for all positions
position = [50 200 850 580];
ms = 8;
ms_dot = 25;
lw = 2;
f1 = figure('Position',position);
Ax1(1) = axes(f1);
glavne_barve = get(Ax1(1), 'ColorOrder');

semilogy(n_homo,t_homo,'.-','Color', glavne_barve(1,:),'MarkerSize',ms_dot,'LineWidth',lw);
hold on
semilogy(n_mep,t_mep,':pentagram','Color', glavne_barve(3,:),'MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw);
hold off
        ylabel('time [s]','FontSize',20)
        xlabel('N','FontSize',20)
        % title('k=3,5,7','FontSize',15)
        ax = gca;
        ax.FontSize=20;
    lgd3 = legend('homotopy','MultiParEig', 'Location', 'southeast','FontSize',20);
axis([4 16 3e-2 4e3])
title('ARMA(2,1)')

