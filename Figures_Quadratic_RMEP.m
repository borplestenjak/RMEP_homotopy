% Plots figures for the paper

% Bor Plestenjak 2026

R = readmatrix("Comparison_polynomial_RMEP.xlsx","Sheet","All data");

% same size for all positions
position = [50 200 850 580];

rows_quad = find(R(:,1)==2);
R = R(rows_quad,2:end); % submatrix with data for quadratic RMEPs

%% ==================================================================
% times for fixed k and different n
set_k = [2 3 4];
set_maxn = [20 12 7];
ms = 8;
lw = 2;

f1 = figure('Position',position);
Ax1(1) = axes(f1);
for j = 1:length(set_k)
    k = set_k(j);
    maxn = set_maxn(j);
    ind = find(R(:,1)==k);
    hom_ind = intersect(find(R(ind,10)>0),find(R(ind,2)<=maxn));
    mep_ind = intersect(find(R(ind,13)>0),find(R(ind,2)<=maxn));
    mac_ind = intersect(find(R(ind,16)>0),find(R(ind,2)<=maxn));
    hom_x = R(ind(hom_ind),2);
    mep_x = R(ind(mep_ind),2);
    mac_x = R(ind(mac_ind),2);
    hom_y0 = R(ind(hom_ind),9);
    hom_y = R(ind(hom_ind),10);
    mep_y = R(ind(mep_ind),13);
    mac_y = R(ind(mac_ind),16);
    
    set(gca,'ColorOrderIndex',j)
    semilogy(mep_x,mep_y,':pentagram','MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw);
    hold on
    set(gca,'ColorOrderIndex',j)
    semilogy(mac_x,mac_y,'--o','MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw);
    set(gca,'ColorOrderIndex',j)
    semilogy(hom_x,hom_y,'.-','MarkerSize',20,'LineWidth',lw);
    ylabel('time [s]','FontSize',20)
    xlabel('n','FontSize',20)
    title('d=2','FontSize',20)
    ax = gca;
    ax.FontSize=20;
end
glavne_barve = get(Ax1(1), 'ColorOrder');
Ax1(2) = copyobj(Ax1(1),gcf);
delete(get(Ax1(2),'children'))
% plot helper data, but invisible
hold on
H1 = plot(nan, nan, '-', 'LineWidth',lw, 'Color', glavne_barve(1,:), 'Parent', Ax1(2),'Visible', 'on');
H2 = plot(nan, nan, '-', 'LineWidth',lw, 'Color', glavne_barve(2,:), 'Parent', Ax1(2),'Visible', 'on');
H3 = plot(nan, nan, '-', 'LineWidth',lw, 'Color', glavne_barve(3,:), 'Parent', Ax1(2),'Visible', 'on');
hold off
% make second axes invisible
set(Ax1(2), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
lgd2 = legend([H1 H2 H3],'k=2', 'k=3', 'k=4', 'Location', 'southeast','FontSize',20);
set(lgd2,'color','none')
set(lgd2, 'TextColor', 'black');
hold off
Ax1(3) = copyobj(Ax1(2),gcf);
delete(get(Ax1(3),'children'))
% plot helper data, but invisible
hold on
G1 = plot(nan, nan, 'k.-','MarkerSize',20,'LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
G2 = plot(nan, nan, 'k--o','MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
G3 = plot(nan, nan, 'k:pentagram','MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
hold off
% make second axes invisible
set(Ax1(3), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
lgd3 = legend([G1 G2 G3], 'homotopy','MacaulayLab','MultiParEig', 'Location', 'northeast','FontSize',20);
set(lgd3,'color','none')
set(lgd3, 'TextColor', 'black');
hold off
axis([1.5 20.5 1.5e-3 3e3])

