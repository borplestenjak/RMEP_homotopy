% Plots figures for the paper

% Bor Plestenjak 2026

R = readmatrix("Comparison_linear_RMEP.xlsx","Sheet","All data");

% same size for all positions
position = [50 200 850 580];

% select just one of the plots 
plotfig = 4; % select figure (1-4)

%% ==================================================================
% first plot 
% times for fixed k and different n
if plotfig==1 
    set_k = [2 3 5 7];
    set_maxn = [50 45 17 12];
    ms = 8;
    ms_dot = 25;
    lw = 2;
    f1 = figure('Position',position);
    Ax1(1) = axes(f1);
    glavne_barve = get(Ax1(1), 'ColorOrder');
    for j = length(set_k):-1:1
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
        semilogy(mac_x,mac_y,'--o','Color', glavne_barve(j,:),'MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw);
        hold on
        semilogy(hom_x,hom_y,'.-','Color', glavne_barve(j,:),'MarkerSize',ms_dot,'LineWidth',lw);
        semilogy(mep_x,mep_y,':pentagram','Color', glavne_barve(j,:),'MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw);
        hold on
        ylabel('time [s]','FontSize',20)
        xlabel('n','FontSize',20)
        ax = gca;
        ax.FontSize=20;
        axis tight
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
    hold off
    % make second axes invisible
    set(Ax1(2), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    lgd2 = legend([H4 H3 H2 H1],'k=7', 'k=5', 'k=3', 'k=2', 'Location', 'southeast','FontSize',20);
    set(lgd2,'color','none')
    set(lgd2, 'TextColor', 'black');
    hold off
    Ax1(3) = copyobj(Ax1(2),gcf);
    delete(get(Ax1(3),'children'))
    % plot helper data, but invisible
    hold on
    G1 = plot(nan, nan, 'k.-','MarkerSize',ms_dot,'LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
    G2 = plot(nan, nan, 'k--o','MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
    G3 = plot(nan, nan, 'k:pentagram','MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
    hold off
    % make second axes invisible
    set(Ax1(3), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    lgd3 = legend([G2 G1 G3], 'MacaulayLab','homotopy','MultiParEig', 'Location', 'south','FontSize',20);
    set(lgd3,'color','none')
    set(lgd3, 'TextColor', 'black');
    hold off
    axis([1 51 3e-4 5e3])
end

%% ==================================================================
% second plot 
% times for fixed n and different k
if plotfig==2 
    set_n = [4 6 8];
    set_maxk = [10 10 8];
    ms = 8;
    ms_dot = 25;
    lw = 2;
    
    f1 = figure('Position',position);
    Ax1(1) = axes(f1);
    glavne_barve = get(Ax1(1), 'ColorOrder');
    for j = 1:length(set_n)
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
        
        semilogy(hom_x,hom_y,'.-','Color', glavne_barve(j,:),'MarkerSize',ms_dot,'LineWidth',lw);
        hold on
        semilogy(mac_x,mac_y,'--o','Color', glavne_barve(j,:),'MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw);
        semilogy(mep_x,mep_y,':pentagram','Color', glavne_barve(j,:),'MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw);
        hold on
        ylabel('time [s]','FontSize',20)
        xlabel('k','FontSize',20)
        ax = gca;
        ax.FontSize=20;
        % axis tight
    end
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
    lgd2 = legend([H3 H2 H1],'n=8', 'n=6', 'n=4', 'Location', 'southeast','FontSize',20);
    set(lgd2,'color','none')
    set(lgd2, 'TextColor', 'black');
    hold off
    Ax1(3) = copyobj(Ax1(2),gcf);
    delete(get(Ax1(3),'children'))
    % plot helper data, but invisible
    hold on
    G1 = plot(nan, nan, 'k.-','MarkerSize',ms_dot,'LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
    G2 = plot(nan, nan, 'k--o','MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
    G3 = plot(nan, nan, 'k:pentagram','MarkerSize',ms,'MarkerFaceColor','w','LineWidth',lw, 'Parent', Ax1(3),'Visible', 'on');
    hold off
    % make second axes invisible
    set(Ax1(3), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    lgd3 = legend([G2 G1 G3], 'MacaulayLab','homotopy','MultiParEig', 'Location', 'south','FontSize',20);
    set(lgd3,'color','none')
    set(lgd3, 'TextColor', 'black');
    hold off
    axis([1.8 10.2 5e-4 1e4])
end

%% ==================================================================
% third plot 
% maximal residual for fixed k and different n
if plotfig==3
    set_k = [2 3];
    set_maxn = [30 30];
    ms = 7;
    lw = 1.5;
    f1 = figure('Position',position);
    Ax1(1) = axes(f1);
    glavne_barve = get(Ax1(1), 'ColorOrder');
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
        hom_y = R(ind(hom_ind),11);
        mep_y = R(ind(mep_ind),14);
        mac_y = R(ind(mac_ind),17);
        
        if j ==1
            p1 = semilogy(hom_x,hom_y,'square','Color', glavne_barve(1,:),'MarkerSize',ms+1,'LineWidth',lw);
        else
            p2 = semilogy(hom_x,hom_y,'square','Color', glavne_barve(1,:),'MarkerSize',ms+1,'MarkerFaceColor','b','LineWidth',lw);
        end
        hold on
        if j ==1
            p3 = semilogy(mep_x,mep_y,'o','Color', glavne_barve(2,:),'MarkerSize',ms,'LineWidth',lw);
            p5 = semilogy(mac_x,mac_y,'d','Color', glavne_barve(3,:),'MarkerSize',ms,'LineWidth',lw);
        else
            p4 = semilogy(mep_x,mep_y,'o','Color', glavne_barve(2,:),'MarkerSize',ms,'MarkerFaceColor',glavne_barve(2,:),'LineWidth',lw);
            p6 = semilogy(mac_x,mac_y,'d','Color', glavne_barve(3,:),'MarkerSize',ms,'MarkerFaceColor',glavne_barve(3,:),'LineWidth',lw);
        end
    end
        legend([p1 p2 p5 p6 p3 p4],'homotopy (k=2)','homotopy (k=3)','MacaulayLab (k=2)','MacaulayLab (k=3)','MultiParEig (k=2)','MultiParEig (k=3)','Location','northwest','FontSize',20)
        ylabel('maximal residual','FontSize',20)
        xlabel('n','FontSize',15)
        ax = gca;
        ax.FontSize=20;
        axis tight
    ax = gca;
    ax.FontSize=20;
    axis([1 31 2e-16 3e-4])
    hold off
end

%% ==================================================================
% fourth plot 
% average steps of homotopy for k=2,3,5
if plotfig==4
    barve = [ 0 0 0; 0 0 0; 6.6000e-02   4.4300e-01   7.4500e-01
      8.6600e-01   3.2900e-01            0
      9.2900e-01   6.9400e-01   1.2500e-01
      5.2100e-01   8.6000e-02   8.1900e-01
      2.3100e-01   6.6600e-01   1.9600e-01
      1.8400e-01   7.4500e-01   9.3700e-01
      8.1900e-01   1.5000e-02   5.4500e-01 ];
    f1 = figure('Position',position);
    set(f1,'defaultAxesColorOrder',barve)
    Ax1(1) = axes(f1);
    glavne_barve = barve(3:end,:);
    yyaxis left
    ms = 7;
    ms_dot = 25;
    lw = 2;
    maxn = 16;
    ind1 = find(R(:,1)==2);
    ind2 = find(R(:,1)==3);
    ind3 = find(R(:,1)==5);
    ind4 = find(R(:,1)==7);
    hom_ind1 = intersect(find(R(ind1,10)>0),find(R(ind1,2)<=maxn));
    hom_ind2 = intersect(find(R(ind2,10)>0),find(R(ind2,2)<=maxn));
    hom_ind3 = intersect(find(R(ind3,10)>0),find(R(ind3,2)<=maxn));
    hom_ind4 = intersect(find(R(ind4,10)>0),find(R(ind4,2)<=maxn));
    hom_x1 = R(ind1(hom_ind1),2);
    hom_x2 = R(ind2(hom_ind2),2);
    hom_x3 = R(ind3(hom_ind3),2);
    hom_x4 = R(ind4(hom_ind4),2);
    hom_y1 = R(ind1(hom_ind1),8);
    hom_y2 = R(ind2(hom_ind2),8);
    hom_y3 = R(ind3(hom_ind3),8);
    hom_y4 = R(ind4(hom_ind4),8);
    
    p1 = plot(hom_x1,hom_y1,'.-','Color', glavne_barve(1,:),'MarkerSize',ms_dot,'LineWidth',lw);
    hold on
    p2 = plot(hom_x2,hom_y2,'.-','Color', glavne_barve(2,:),'MarkerSize',ms_dot,'LineWidth',lw);
    p3 = plot(hom_x3,hom_y3,'.-','Color', glavne_barve(3,:),'MarkerSize',ms_dot,'LineWidth',lw);
    p4 = plot(hom_x4,hom_y4,'.-','Color', glavne_barve(4,:),'MarkerSize',ms_dot,'LineWidth',lw);
    ylabel('average steps','FontSize',20)
    xlabel('n','FontSize',20)
    ax = gca;
    ax.FontSize=20;    
    % axis tight
    axis([1 17 -2 60])
    xticks(2:2:16)
    yyaxis right
    ind1 = find(R(:,1)==2);
    ind2 = find(R(:,1)==3);
    ind3 = find(R(:,1)==5);
    hom_x1 = R(ind1(hom_ind1),2);
    hom_x2 = R(ind2(hom_ind2),2);
    hom_x3 = R(ind3(hom_ind3),2);
    hom_x4 = R(ind4(hom_ind4),2);
    hom_y1 = 100*R(ind1(hom_ind1),7);
    hom_y2 = 100*R(ind2(hom_ind2),7);
    hom_y3 = 100*R(ind3(hom_ind3),7);
    hom_y4 = 100*R(ind4(hom_ind4),7);
    
    r1 = plot(hom_x1,hom_y1,'square:','Color', glavne_barve(1,:),'MarkerSize',ms+1,'LineWidth',lw);
    hold on
    r2 = plot(hom_x2,hom_y2,'square:','Color', glavne_barve(2,:),'MarkerSize',ms+1,'LineWidth',lw);
    r3 = plot(hom_x3,hom_y3,'square:','Color', glavne_barve(3,:),'MarkerSize',ms+1,'LineWidth',lw);
    r4 = plot(hom_x4,hom_y4,'square:','Color', glavne_barve(4,:),'MarkerSize',ms+1,'LineWidth',lw);
    ylabel('retraced paths [%]','FontSize',20)

    legend([p1 p2 p3 p4 r1 r2 r3 r4],'k=2: steps','k=3: steps','k=5: steps','k=7: steps','retr. paths','retr. paths','retr. paths','retr. paths','Location','northwest','FontSize',20,'NumColumns',2)
    axis([1 17 -1 25])
    hold off
end

