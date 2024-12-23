%Create Figure 1 DGH Trajectories
%Assuming activating from the main folder
current_cd=cd;
load(horzcat(current_cd,'\Data\Fig1DGH\WS_Fig1.mat'));
%This are the trajectories taken from -2.5 K_BT simulation that also
%produced the movies


figgg=figure;
yyaxis left
bbb=get(figgg,'Position');
h_factor=bbb(3)/bbb(4);
new_width=8.7;
% load('Built_mat_discretePLoopaUP_1_time_0_85_amp.mat');
% Median_mat_discrete(2)=Median_mat_discrete(2)/10;
set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 4/3*new_width]);

x=New_Times;
subplot(4,1,1);
%imshow(I, 'border', 'tight');
energy(1)=-7.5;
plot(x,energy);
xlim([0 x(end)]);
hold on;
scatter(x(1), energy(1), 10, 'filled','MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
set(gca,'Xticklabel',[])
% set(gca,'FontSize',24)
% y=ylabel(horzcat('$E \ ', ' [k_{B}T]$'),'FontSize',6,'Interpreter','latex');
y= ylabel( '{\boldmath${E\ [k_{B}T]}$}','Interpreter','latex');
set(gca,'FontSize',6);

subplot(4,1,3);
% set(y, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% axes3=axes(f, 'Position', positionVector3);
%imshow(I, 'border', 'tight');
plot(x,dist1);
xlim([0 x(end)]);
set(gca,'Xticklabel',[])
y3=ylabel('{\boldmath${d_{1}}$}','Interpreter','latex');
% set(y3, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% % y3.Position(1)=0.623958333333333;
% set(gca,'FontSize',24)
set(gca,'FontSize',6);

hold on;
subplot(4,1,4);
% axes4=axes(f, 'Position', positionVector4);
%imshow(I, 'border', 'tight');
plot(x,dist2);
xlim([0 x(end)]);
% xlabel({'Time[s]',' '},'FontSize',6)
% xx0= xlabel( '{},{\boldmath${Time\ [s]}$}','Interpreter','latex');
xx0= xlabel({'$$\textbf {\hspace{0.1mm} }$$', '$$\textbf {Time [s]}$$'},'Interpreter','latex');
% y4=ylabel('$d_{2}$  ','FontSize',6,'Interpreter','latex');
y4=ylabel('{\boldmath${d_{2}}$}','Interpreter','latex');
% set(y4, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
set(gca,'FontSize',6);
% cd('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC');
nameXTotal='Fig1D';
print (nameXTotal,'-dpng','-r600');
% cd('C:\Users\admin\Pictures');
% print (nameXTotal,'-dpng','-r600');

