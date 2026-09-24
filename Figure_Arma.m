n_homo = [4:2:34]';

t_homo = [
 0.739 
 1.352 
 2.995 
 6.300 
 11.462 
 21.009 
 35.285 
 55.040 
 72.686 
 108.015 
 149.546 
 205.617 
 276.834 
 367.879 
 569.7886
 776.5583
];

n_mep = [4:2:34]';

t_mep = [
 0.006 
 0.014 
 0.030 
 0.062 
 0.152 
 0.275 
 0.715 
 1.480 
 2.890 
 5.400 
 9.140 
 15.134 
 24.049 
 37.215 
 53.883 
 78.073 
 ];

n_mac = [4:2:10]';

t_mac = [
 0.409 
 7.573 
98.502
1204.8
];

% same size for all positions
position = [50 200 850 580];
ms = 8;
ms_dot = 25;
lw = 2;
f1 = figure('Position',position);
Ax1(1) = axes(f1);
glavne_barve = get(Ax1(1), 'ColorOrder');

semilogy(n_mac,t_mac,'--o','Color', glavne_barve(2,:),'MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw);
hold on
semilogy(n_homo,t_homo,'.-','Color', glavne_barve(1,:),'MarkerSize',ms_dot,'LineWidth',lw);
semilogy(n_mep,t_mep,':pentagram','Color', glavne_barve(3,:),'MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw);
hold off
        ylabel('time [s]','FontSize',20)
        xlabel('N','FontSize',20)
        % title('k=3,5,7','FontSize',15)
        ax = gca;
        ax.FontSize=20;
    lgd3 = legend('MacaulayLab','homotopy','MultiParEig', 'Location', 'southeast','FontSize',20);
axis([3 36 3e-3 3e3])
title('ARMA(1,1)')
