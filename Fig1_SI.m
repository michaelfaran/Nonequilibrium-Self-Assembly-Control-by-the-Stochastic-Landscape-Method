 % % load('WS1.mat');
 clear;
 % cd('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC');
 % clear;
 %  addpath ('C:\Users\admin\Documents\GitHub\GradProject');
 current_cd=cd;
addpath(horzcat(cd,'\Data\SI_Fig1_SI_Fig2'));
energies_true=[2 2.2 2.4 2.5 2.6 2.8 3  3.2  3.4 3.5 3.6];% 3.2 3.4
energies_true=[2 2.2 2.4 2.5 2.6 2.8 3 3.1 3.2 3.3 3.4 3.5 3.6];
Aa=[ 3 4];
Aa2=[ 25 25];
figgg=figure;
bbb=get(figgg,'Position');
h_factor=bbb(3)/bbb(4);
new_width=8.7;
% load('Built_mat_discretePLoopaUP_1_time_0_85_amp.mat');
% Median_mat_discrete(2)=Median_mat_discrete(2)/10;
set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width (length(Aa)/2*new_width)*1.1]);
helper=length(Aa2);
for kk=1:1:helper
subplot(helper,1,kk);
TN=horzcat(num2str(Aa(kk)),'_T_',num2str(Aa2(kk)),'_N');
Adressme=horzcat('mimi_SS_god_',TN,'_eqUP.mat');
load(Adressme);
%try
%Adressme2=horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\Built_mat_discretePLoopaUP_Eq_',regexprep(TN,'N','P'),'.mat');
%load(Adressme2);
%catch
Adressme2=horzcat('Built_mat_discretePLoopaUP_N_Control_C_',regexprep(TN,'N','P'),'.mat');
load(Adressme2);
%end
[X, Y] = extract_X_Y_N(Adressme);
load(horzcat('mimi_SS_god_discrete_',num2str(X),'_T_',num2str(Y),'_N_eqUP.mat'));
J=energies_true;
ZZ=squeeze(mimi_SS_god);
ZZ_where_eq=squeeze(mimi_SS_god_discrete);
% ZZ(2,:)=ZZ(2,:)/10;
maxZ_vec=zeros(1,length(J));
for yy=1:1:(length(J)) 
maxZ=max(ZZ(yy,:));
% maxZ_vec=maxZ;
ZZ(yy,find(ZZ_where_eq(yy,:)==2*10^7))=maxZ;
ZZ(yy,:)=ZZ(yy,:)./maxZ;
end

Mmm=median(ZZ,2);
% Mmm(2)=Mmm(2)/10;

% figgg=figure;
yyaxis left
% bbb=get(figgg,'Position');
% h_factor=bbb(3)/bbb(4);
% new_width=8.7;
% load('Built_mat_discretePLoopaUP_1_time_0_85_amp.mat');
% Median_mat_discrete(2)=Median_mat_discrete(2)/10;
% set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);
%energies_true=[2 2.2 2.4 2.5 2.6 2.8 3 3.1 3.2 3.3 3.4 3.5];% 3.2 3.4
% plot(flip(-energies_true,2),flip(Median_mat_discrete,1),'.b');
% hold on;
% plot(flip(-energies_true,2),flip(Mmm,1),'.r');
% ylabel('Discrete Med(Tfas)')
% xlabel('Js[k_{B}T]')
% print (figgg,'-dpng','-r300');
% 
% aa=subplot(2,1,1);
pop=Mmm;
% figure; scatter (J,flip(ZZ),60,'filled');
ylabel('{\boldmath${\hat{T}_{FAS}}$}','FontSize',6,'interpreter','latex')
J=energies_true;
% xlabel('$J_{s}$','FontSize',6,'interpreter','latex')
hold on;
plot(-flip(J,2),flip(pop',2),'Color',  [0.5 0.5 1],'LineWidth',2,'LineStyle','--');
hold on;
% DOGG='Eq vs Drive Plot, continous';
% 
% ZZ=squeeze(mimi_SS_god);
% ZZ_where_eq=squeeze(mimi_SS_god_discrete);
% % ZZ(2,:)=ZZ(2,:)/10;
% maxZ_vec=zeros(1,length(J));
% for yy=1:1:(length(J)) 
% maxZ=max(ZZ(yy,:));
% maxZ_vec=maxZ;
% ZZ(yy,find(ZZ_where_eq(yy,:)==2*10^7))=maxZ;
% end
% 
% ZZ=ZZ./maxZ;

for yy=1:1:(length(J)) 
    a=(ZZ(yy,:));
[C,ia,ic] = unique(a);
a_counts = accumarray(ic,1);
try
value_counts = [C, a_counts];
catch
value_counts = [C, a_counts'];

end
% BubblePlot(ones(1,size(a_counts,1))*J(length(J)-yy+1), C,0.5*log(a_counts)+1,'b',8);
dd=colormap(cool);%following gili request for xomplemntary
bubbleplotUP(-ones(1,size(a_counts,1))*J(yy), C,0.5*log(a_counts)+1,dd(1,:),8,1,sum(a_counts),max(a_counts),size(ZZ,2));
hold on;
% BubblePlot(ones(1,size(1,1))*J(yy), C(end),0.5*log(a_counts(end))+1,'k',8);
% hold on;
end
if kk==1
cb = colorbar('northoutside');
cbscalein = cb.Limits;
cbscaleout = [0 100];
cb.Ticks=0:0.2:1.1;
cb.TickLabels = [{'0%','20%','40%','60%','80%','100%'}];
end
% title('{\boldmath${\hat{T}_{FAS}(Js),\, Eq}$}','FontSize',6,'interpreter','latex')
 % ylim([0 105]);
set(gca,'FontSize',6)
hold on;
% bb=subplot(2,1,2);
% 
% cd(Rafiq_of_the_many);
% Heated_debate_map_CALL(Main_dir_hazirim,energies_true,drive_vec,Built_mat_discrete,Median_mat_discrete,Median_mat);
% Built_mat_discreteP=Built_mat_discrete./turncoat;
% path_str =cd;
% cd('C:\Users\admin\Documents\RMFO2');
% % Split the string by backslashes
% parts = strsplit(path_str, '\');
% 
% % Get the last part
% last_part = parts{end};
% 
% save(horzcat('Built_mat_discreteP',last_part),"Built_mat_discreteP")
% cd(Rafiq_of_the_many);
% addpath 'C:\Users\admin\Documents\RMFO2'
% figgg=figure;
% bbb=get(figgg,'Position');
% h_factor=bbb(3)/bbb(4);
% new_width=8.7;
% load('Built_mat_discretePLoopaUP_N_Control_New_.mat');
yyaxis right
% set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);
%energies_true=[2 2.2 2.4 2.5 2.6 2.8 3 3.1 3.2 3.3 3.4 3.5];% 3.2 3.4
plot(flip(-energies_true,2),100*flip(Built_mat_discreteP,1),'ok-','MarkerSize',4);
% load('Built_mat_discretePLoopaUP_1_time_0_85_amp2.mat');
% hold on;
% plot(flip(-energies_true,2),100*flip(Built_mat_discreteP,1),'.r--','MarkerSize',6,'Marker','square');
ylabel('{\boldmath${SA[\%]}$}','FontSize',6,'interpreter','latex');
title(horzcat('{\boldmath${M_{T}=',num2str(Aa(kk)),',N=',num2str(Aa2(kk)),'}$}'),'FontSize',6,'interpreter','latex');
if kk==helper
xlabel('{\boldmath${J_{s}[k_{B}T]}$}','FontSize',6,'interpreter','latex');
end
ylim([-10 105]);
set(gca,'YTick',0:20:100);

ax=gca;
ax.YAxis(1).Color = [0.5 0.5 1];
ax.YAxis(2).Color = 'k';

% gca.Ticks=0:20:100;
% gca.TickLabels = [{'0%','20%','40%','60%','80%','100%'}];
% h=legend('$Eq$','$\Delta\rho=1.5,W=9\cdot10^4$','interpreter','latex');
% h.ItemTokenSize(1)=20;
% set(h,...
%     'Position',[0.333668902747363 0.162799085309368 0.31808942362545 0.216463409000781]);

% DOGG=horzcat('KMC_TFAS_vs_Js_SA_w_drive_',TN);
% ax=gca;
%  ax.YAxis(1).Color = [0.5 0.5 1];
% ax.YAxis(2).Color = 'k';
% set(gca,'Position',[0.1300 0.1716 0.7750 0.5426]);
% % % 
%  print (DOGG,'-dpng','-r300');
% cd(Rafiq_of_the_many);
hold on;
end
Y=0.0182*2/length(Aa2);
cb.Position=[0.1311    0.9291-Y   0.7744    Y];
cb.Units='pixels';
pause( 1);
XXX=cb.Position;
Y=XXX(4);
% cb.Position=[0.1311    0.9291    0.7744    Y];

% Find all the axes in the current figure
fig = gcf;  % Get the current figure handle
all_axes = findall(fig, 'Type', 'axes');  % Find all axes in the figure
adjustment = [0, 0, 0, Y];

% Loop through each axis and adjust the position
for ii=1:1:length(all_axes)
ax = all_axes(ii);  % Get the current axis
% current_pos = get(ax, 'Position');  % Get current axis position
set(ax,'Units','pixels');
pause( 1);
current_pos = get(ax, 'Position'); 
pause( 1);
if ii~=length(all_axes)
 adjustment = [0, -Y, 0, Y];
else
 adjustment = [0, 0, 0, Y];
   
end
CC=current_pos - 4*adjustment;
%You have toa ctivate ehre if you want figure to be clear.. stupid issue
set(ax, 'Position', CC);
pause( 1);
% Subtract the adjustment
end
% set(gca,'FontSize',6);
if sum(Aa2)==length(Aa2)*25
DOGG='Fig2EXTSI';
text(-0.07,2.57, 'A', 'Units', 'normalized', 'FontSize', 18, 'FontName', 'Calibri', 'HorizontalAlignment', 'center');
text(-0.07, 1.1, 'B', 'Units', 'normalized', 'FontSize', 18, 'FontName', 'Calibri', 'HorizontalAlignment', 'center');
else
DOGG='Fig2EXTSI_36';
end
% cd('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC');
print (DOGG,'-dpng','-r300');
% cd('C:\Users\admin\Pictures');
% print (DOGG,'-dpng','-r300');