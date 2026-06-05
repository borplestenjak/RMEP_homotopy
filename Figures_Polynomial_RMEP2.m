% Plots figures for the paper

% Bor Plestenjak 2026

Rall = readmatrix("Comparison_polynomial_RMEP.xlsx","Sheet","All data");

% same size for all positions
position = [50 200 850 580];
plot1 = 1;
plot2 = 0;
plot3 = 0;
plot4 = 0;
plot5 = 0;

marker_values = [2, 3, 4, 5, 6, 7];
unicode_decimals = 9311 + marker_values; 
marker_strings = cellstr(char(unicode_decimals.'));

rows_3 = find(Rall(:,1)==7);
R = Rall(rows_3,2:end); % submatrix with data for cubic RMEPs

%% ==================================================================
% first plot 
% times for fixed k and different n
if plot1 
    set_deg = [3 4 5 6 7];
    set_maxn = [10 10 10 10 10];
c    for j = 1:length(set_deg)
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
        semilogy(hom_x,hom_y,'.-','MarkerSize',ms,'LineWidth',lw);
        hold on
        set(gca,'ColorOrderIndex',j)
        semilogy(mac_x,mac_y,'--o','MarkerSize',6,'MarkerFaceColor','w','LineWidth',lw);
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
    lgd2 = legend([H1 H2 H3 H4 H5], 'd=3', 'd=4', 'd=5', 'd=6', 'd=7', 'Location', 'southeast','FontSize',15);
    set(lgd2,'color','none')
    set(lgd2, 'TextColor', 'black');
    hold off
    Ax1(3) = copyobj(Ax1(2),gcf);
    delete(get(Ax1(3),'children'))
    % plot helper data, but invisible
    hold on
    G1 = plot(nan, nan, 'k.-','MarkerSize',ms,'LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
    G2 = plot(nan, nan, 'k--o','MarkerSize',6,'MarkerFaceColor','w','LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
    hold off
    % make second axes invisible
    set(Ax1(3), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    lgd3 = legend([G1 G2], 'homotopy','Macaulay', 'Location', 'northwest','FontSize',15);
    set(lgd3,'color','none')
    set(lgd3, 'TextColor', 'black');
    hold off
end
    
%% ==================================================================
% second plot 
% times for fixed n and different k
if plot2 
    set_n = [2 4 6];
    set_maxk = [4 4 4];
    
    for j = 1:length(set_k)
        figure('Position',position)
        n = set_n(j);
        maxk = set_maxk(j);
        ind = find(R(:,2)==n);
        hom_ind = intersect(find(R(ind,10)>0),find(R(ind,1)<=maxk));
        mep_ind = intersect(find(R(ind,13)>0),find(R(ind,1)<=maxk));
        mac_ind = intersect(find(R(ind,16)>0),find(R(ind,1)<=maxk));
        hom_x = R(ind(hom_ind),1);
        mep_x = R(ind(mep_ind),1);
        mac_x = R(ind(mac_ind),1);
        hom_y0 = R(ind(hom_ind),9);
        hom_y = R(ind(hom_ind),10);
        mep_y = R(ind(mep_ind),13);
        mac_y = R(ind(mac_ind),16);
        
        semilogy(hom_x,hom_y,'b-square','MarkerSize',10,'MarkerFaceColor','b','LineWidth',3);
        hold on
        semilogy(mep_x,mep_y,'r:o','MarkerSize',10,'MarkerFaceColor','r','LineWidth',3);
        semilogy(mac_x,mac_y,'g--d','MarkerSize',10,'MarkerFaceColor','g','LineWidth',3 );
        semilogy(hom_x,hom_y0,'b square','MarkerSize',10,'LineWidth',2);
        hold off
        legend('homotopy','MultiParEig','Macaulay','Location','southeast','FontSize',15)
        ylabel('time [s]','FontSize',15)
        xlabel('n','FontSize',15)
        title(sprintf('n=%d',n),'FontSize',15)
        ax = gca;
        ax.FontSize=20;
        axis tight
    end
end

%% ==================================================================
% third plot 
% maximal residual for fixed k and different n
if plot3
    set_k = [2 3];
    set_maxn = [20 20];
    
    for j = 1:length(set_k)
        figure('Position',position)
        k = set_k(j);
        maxn = set_maxn(j);
        ind = find(R(:,1)==k);
        hom_ind = intersect(find(R(ind,10)>0),find(R(ind,2)<=maxn));
        mep_ind = intersect(find(R(ind,13)>0),find(R(ind,2)<=maxn));
        mac_ind = intersect(find(R(ind,16)>0),find(R(ind,2)<=maxn));
        hom_x = R(ind(hom_ind),2);
        mep_x = R(ind(mep_ind),2);
        mac_x = R(ind(mac_ind),2);
        hom_y = R(ind(hom_ind),11);
        mep_y = R(ind(mep_ind),14);
        mac_y = R(ind(mac_ind),17);
        
        semilogy(hom_x,hom_y,'b square','MarkerSize',10,'LineWidth',2);
        hold on
        semilogy(mep_x,mep_y,'r o','MarkerSize',10,'LineWidth',2);
        semilogy(mac_x,mac_y,'g d','MarkerSize',10,'LineWidth',2);
        hold off
        legend('homotopy','MultiParEig','Macaulay','Location','northwest','FontSize',15)
        ylabel('maximal residual','FontSize',15)
        xlabel('n','FontSize',15)
        title(sprintf('k=%d',k),'FontSize',15)
        ax = gca;
        ax.FontSize=20;
        axis tight
    end
end

%% ==================================================================
% fourth plot 
% average steps of homotopy for k=2,3,5
if plot4
    maxn = 15;
    ind1 = find(R(:,1)==2);
    ind2 = find(R(:,1)==3);
    ind3 = find(R(:,1)==5);
    hom_ind1 = intersect(find(R(ind1,10)>0),find(R(ind1,2)<=maxn));
    hom_ind2 = intersect(find(R(ind2,10)>0),find(R(ind2,2)<=maxn));
    hom_ind3 = intersect(find(R(ind3,10)>0),find(R(ind3,2)<=maxn));
    hom_x1 = R(ind1(hom_ind1),2);
    hom_x2 = R(ind2(hom_ind2),2);
    hom_x3 = R(ind3(hom_ind3),2);
    hom_y1 = R(ind1(hom_ind1),8);
    hom_y2 = R(ind2(hom_ind2),8);
    hom_y3 = R(ind3(hom_ind3),8);
    
    plot(hom_x1,hom_y1,'b square','MarkerSize',10,'LineWidth',2);
    hold on
    plot(hom_x2,hom_y2,'r o','MarkerSize',10,'LineWidth',2);
    plot(hom_x3,hom_y3,'g d','MarkerSize',10,'LineWidth',2 );
    hold off
    legend('k=2','k=3','k=5','Location','northwest','FontSize',15)
    ylabel('average steps','FontSize',15)
    xlabel('n','FontSize',15)
    ax = gca;
    ax.FontSize=15;    
    axis tight
    xticks(2:2:16)
end

%% ==================================================================
% fifth plot 
% percentage of repeated steps of homotopy for k=2,3,5
if plot5
    maxn = 15;
    ind1 = find(R(:,1)==2);
    ind2 = find(R(:,1)==3);
    ind3 = find(R(:,1)==5);
    hom_ind1 = intersect(find(R(ind1,10)>0),find(R(ind1,2)<=maxn));
    hom_ind2 = intersect(find(R(ind2,10)>0),find(R(ind2,2)<=maxn));
    hom_ind3 = intersect(find(R(ind3,10)>0),find(R(ind3,2)<=maxn));
    hom_x1 = R(ind1(hom_ind1),2);
    hom_x2 = R(ind2(hom_ind2),2);
    hom_x3 = R(ind3(hom_ind3),2);
    hom_y1 = 100*R(ind1(hom_ind1),7);
    hom_y2 = 100*R(ind2(hom_ind2),7);
    hom_y3 = 100*R(ind3(hom_ind3),7);
    
    plot(hom_x1,hom_y1,'b square','MarkerSize',10,'LineWidth',2);
    hold on
    plot(hom_x2,hom_y2,'r o','MarkerSize',10,'LineWidth',2);
    plot(hom_x3,hom_y3,'g d','MarkerSize',10,'LineWidth',2 );
    hold off
    legend('k=2','k=3','k=5','Location','northwest','FontSize',15)
    ylabel('repeated paths [%]','FontSize',15)
    xlabel('n','FontSize',15)
    ax = gca;
    ax.FontSize=15;
    axis tight
    xticks(2:2:16)
end
