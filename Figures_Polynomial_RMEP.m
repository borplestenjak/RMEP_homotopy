% Plots figures for the paper

% Bor Plestenjak 2026

Rall = readmatrix("Comparison_polynomial_RMEP.xlsx","Sheet","All data");

% same size for all positions
position = [50 200 850 580];

marker_values = [2, 3, 4, 5, 6, 7];
unicode_decimals = 9311 + marker_values; 
marker_strings = cellstr(char(unicode_decimals.'));

rows_3 = find(Rall(:,1)==7);
R = Rall(rows_3,2:end); % submatrix with data for cubic RMEPs

set_deg = [3 4 5 6 7];
set_maxn = [10 10 10 10 10];
f1 = figure('Position',position);
Ax1(1) = axes(f1);
for j = 1:length(set_deg)
    deg = set_deg(j);
    rows = find(Rall(:,1)==deg);
    R = Rall(rows,2:end); % submatrix with data for cubic RMEPs
    maxn = set_maxn(j);
    ind = find(R(:,1)==2);
    hom_ind = intersect(find(R(ind,10)>0),find(R(ind,2)<=maxn));
    mac_ind = intersect(find(R(ind,16)>0),find(R(ind,2)<=maxn));
    hom_x = R(ind(hom_ind),2);
    mac_x = R(ind(mac_ind),2);
    hom_y0 = R(ind(hom_ind),9);
    hom_y = R(ind(hom_ind),10);
    mac_y = R(ind(mac_ind),16);
    set(gca,'ColorOrderIndex',j)
    semilogy(hom_x,hom_y,'.-','MarkerSize',25,'LineWidth',lw);
    hold on
    set(gca,'ColorOrderIndex',j)
    semilogy(mac_x,mac_y,'--o','MarkerSize',7,'MarkerFaceColor','w','LineWidth',lw);
    % legend('homotopy','Macaulay','Location','northwest','FontSize',15)
    ylabel('time [s]','FontSize',15)
    xlabel('n','FontSize',15)
    title('k=2','FontSize',15)
    ax = gca;
    ax.FontSize=20;
    axis([1.8 10.2 2e-2 3e2])
end
glavne_barve = get(Ax1(1), 'ColorOrder');
Ax1(2) = copyobj(Ax1(1),gcf);
delete(get(Ax1(2),'children'))
% plot helper data, but invisible
hold on
H1 = plot(nan, nan, '-', 'LineWidth',lw, 'Color', glavne_barve(1,:), 'Parent', Ax1(2),'Visible', 'on');
H2 = plot(nan, nan, '-', 'LineWidth',lw, 'Color', glavne_barve(2,:), 'Parent', Ax1(2),'Visible', 'on');
H3 = plot(nan, nan, '-', 'LineWidth',lw, 'Color', glavne_barve(3,:), 'Parent', Ax1(2),'Visible', 'on');
H4 = plot(nan, nan, '-', 'LineWidth',lw, 'Color', glavne_barve(4,:), 'Parent', Ax1(2),'Visible', 'on');
H5 = plot(nan, nan, '-', 'LineWidth',lw, 'Color', glavne_barve(5,:), 'Parent', Ax1(2),'Visible', 'on');
hold off
% make second axes invisible
set(Ax1(2), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
lgd2 = legend([H1 H2 H3 H4 H5], 'd=3', 'd=4', 'd=5', 'd=6', 'd=7', 'Location', 'southeast','FontSize',20);
set(lgd2,'color','none')
set(lgd2, 'TextColor', 'black');
hold off
Ax1(3) = copyobj(Ax1(2),gcf);
delete(get(Ax1(3),'children'))
% plot helper data, but invisible
hold on
G1 = plot(nan, nan, 'k.-','MarkerSize',23,'LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
G2 = plot(nan, nan, 'k--o','MarkerSize',6,'MarkerFaceColor','w','LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
hold off
% make second axes invisible
set(Ax1(3), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
lgd3 = legend([G1 G2], 'homotopy','MacaulayLab', 'Location', 'northwest','FontSize',20);
set(lgd3,'color','none')
set(lgd3, 'TextColor', 'black');
hold off
