% Introductory example with 3 x 2 two-parameter linear RMEP and simple real
% eigenvalues. If we start with a real problem, then real homotopy paths
% get stuck at bifurcation points. This shows why we should use complex
% paths even if the proble is real and has real eigenvalues.

% Bor Plestenjak 2026

A0 = [-1 1; 
      2 -2; 
      3 -1];

A1 = [1 2; 
      -2 3; 
      3 3];

A2 = [1 1;
      2 -2;
      -3 4];

% Multipareig 
fprintf('\nMultiParEig\n-----------\n')
opts = [];
lambda = rect_multipareig({A0,A1,A2},opts)

% MacaulayLab 
fprintf('\nMacaulayLab \n-----------\n')
A = {A0, A1, A2};
suppA = [0 0; 1 0; 0 1];
try
    lambda1 = rect_multipareig_macaulay(A,30,true)
catch ME
    fprintf('Error in macaulaylab: %s \n',ME.message)   
end

% Real homotopy 
fprintf('\nHomotopy from an initial real problem and real paths \n----------\n')
B0A = [2 0; -1 1; 0 2];
B1 =  [3 0; 0 1; 0 0];
B2 =  [0 0; 1 0; 0 1];
BA = {B0A,B1,B2};

[Lambda0,X0] = rect_multipareig(BA)
opts = [];
opts.gamma = 1;
opts.display = 1;
[lambdaT, XT, tcell1, ycell1, stat] = homotopy_rmep(A,BA,Lambda0,X0,opts);
lambdaN1 = lambdaT(:,2:end)./lambdaT(:,1)

fprintf('\nHomotopy from an initial complex problem and complex paths  \n----------\n')
B0B = [2+1i 0; 3-1i 1+1i; 0 2-1i];
B1 =  [1 0; 0 1; 0 0];
B2 =  [0 0; 1 0; 0 1];
BB = {B0B,B1,B2};

[Lambda0,X0] = rect_multipareig(BB)
opts = [];
opts.gamma = sqrt(2)/2*(1+1i);
opts.display = 1;
[lambdaT, XT, tcell2, ycell2, stat] = homotopy_rmep(A,BB,Lambda0,X0,opts);
lambdaN2 = lambdaT(:,2:end)./lambdaT(:,1)


%% Figure with paths from the real homotopy
position = [50 200 850 580];
lw = 3;
f1 = figure('Position',position);
Ax1(1) = axes(f1);
glavne_barve = get(Ax1(1), 'ColorOrder');
for j = 1:3
    plot(tcell1{j},real(ycell1{j}(4,:)./ycell1{j}(3,:)),'Color',glavne_barve(j,:),'LineWidth',lw) 
    hold on
end
for j = 1:3
    plot(tcell1{j},real(ycell1{j}(5,:)./ycell1{j}(3,:)),'--','Color',glavne_barve(j,:),'LineWidth',lw) 
end
hold off
legend('\lambda^{(1)}_1','\lambda^{(2)}_1','\lambda^{(3)}_1','\lambda^{(1)}_2','\lambda^{(2)}_2','\lambda^{(3)}_2',...
    'NumColumns',2,'Location','southeast','FontSize',20)
ax = gca;
ax.FontSize=20;
xlabel('t','FontSize',20)
ylabel('Real(\lambda)','FontSize',fs)

%% Figure with paths from the complex homotopy

f2 = figure('Position',position);
Ax2(1) = axes(f2);
glavne_barve = get(Ax2(1), 'ColorOrder');

j=1; pl1 = plot3(real(ycell2{j}(4,:)./ycell2{j}(3,:)),imag(ycell2{j}(4,:)./ycell2{j}(3,:)),tcell2{j},'Color',glavne_barve(j,:),'LineWidth',lw);
hold on
j=2; pl3 = plot3(real(ycell2{j}(4,:)./ycell2{j}(3,:)),imag(ycell2{j}(4,:)./ycell2{j}(3,:)),tcell2{j},'Color',glavne_barve(j,:),'LineWidth',lw);
j=3; pl5 = plot3(real(ycell2{j}(4,:)./ycell2{j}(3,:)),imag(ycell2{j}(4,:)./ycell2{j}(3,:)),tcell2{j},'Color',glavne_barve(j,:),'LineWidth',lw);
j=1; pl2 = plot3(real(ycell2{j}(5,:)./ycell2{j}(3,:)),imag(ycell2{j}(5,:)./ycell2{j}(3,:)),tcell2{j},'--','Color',glavne_barve(j,:),'LineWidth',lw);
j=2; pl4 = plot3(real(ycell2{j}(5,:)./ycell2{j}(3,:)),imag(ycell2{j}(5,:)./ycell2{j}(3,:)),tcell2{j},'--','Color',glavne_barve(j,:),'LineWidth',lw);
j=3; pl6 = plot3(real(ycell2{j}(5,:)./ycell2{j}(3,:)),imag(ycell2{j}(5,:)./ycell2{j}(3,:)),tcell2{j},'--','Color',glavne_barve(j,:),'LineWidth',lw);
ax = gca;
fs = 20;
ax.FontSize=fs;
xlabel('Real(\lambda)','FontSize',fs)
ylabel('Imag(\lambda)','FontSize',fs)
zlabel('t','FontSize',fs)

% Limits in the complex plane
xr = [-3.5, 3.5];
yr = [-3.5, 3.5];

% Mesh for the two complex planes
[X,Y] = meshgrid(linspace(xr(1),xr(2),25), ...
                 linspace(yr(1),yr(2),25));

Z0 = zeros(size(X));
Z1 = ones(size(X));

% Planes t = 0 and t = 1
sr1 = surf(X,Y,Z0, ...
    'FaceAlpha',0.12, ...
    'EdgeAlpha',0.25);

sr2 = surf(X,Y,Z1, ...
    'FaceAlpha',0.12, ...
    'EdgeAlpha',0.25);

for j = 1:3
    plot3(real(ycell2{j}(4,1)./ycell2{j}(3,1)),imag(ycell2{j}(4,1)./ycell2{j}(3,1)),0,'.','Color',glavne_barve(j,:),'MarkerSize',30)
    plot3(real(ycell2{j}(5,1)./ycell2{j}(3,1)),imag(ycell2{j}(5,1)./ycell2{j}(3,1)),0,'.','Color',glavne_barve(j,:),'MarkerSize',30)
end
for j = 1:3
    plot3(real(ycell2{j}(4,end)./ycell2{j}(3,end)),imag(ycell2{j}(4,end)./ycell2{j}(3,end)),1,'.','Color',glavne_barve(j,:),'MarkerSize',30)
    plot3(real(ycell2{j}(5,end)./ycell2{j}(3,end)),imag(ycell2{j}(5,end)./ycell2{j}(3,end)),1,'.','Color',glavne_barve(j,:),'MarkerSize',30)
end

legend([pl1 pl2 pl3 pl4 pl5 pl6 sr2 sr1],{'\lambda^{(1)}_1','\lambda^{(2)}_1','\lambda^{(3)}_1','\lambda^{(1)}_2','\lambda^{(2)}_2','\lambda^{(3)}_2','t=1','t=0'},...
    'NumColumns',1,'Location','northwest','FontSize',fs)

axis([-3.5 3.5 -3.5 3.5 0 1])