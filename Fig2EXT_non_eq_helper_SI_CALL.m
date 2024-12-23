function [EQ_Built_mat_discreteP2,Built_mat_discreteP,Mmm,Mmm_drive]=Fig2EXT_non_eq_helper_SI_CALL(Aa,Aa2,energies_true,Rafiq_of_the_many,helper_length,energies_val)
% load('WS1.mat');
 % clear;
%   addpath ('C:\Users\admin\Documents\GitHub\GradProject');
% energies_true=[2 2.2 2.4 2.5 2.6 2.8 3  3.2  3.4 3.5];% 3.2 3.4
% energies_true=[2 2.2 2.4 2.5 2.6 2.8 3 3.1 3.2 3.3 3.4 3.5];
% 
% Rafiq_of_the_many1='C:\Users\admin\Documents\RMFO\LoopaUP_Aq_Nw_1_t_0_85_A_3_T_25_P';
% Rafiq_of_the_many2='C:\Users\admin\Documents\RMFO\LoopaUP_Aq_Nw_1_t_0_85_A_2_T_36_P';
% Rafiq_of_the_many1='C:\Users\admin\Documents\RMFO\LoopaUP_Ama_1_t_0_85_A_4_T_25_P';
% Rafiq_of_the_many={Rafiq_of_the_many1};% Rafiq_of_the_many2
% energies_true=[ 3.4 3.5 3.6];
% 
% 
% Aa=[3 2];
% Aa2=[25 36];
% 
% Aa=[4];
% Aa2=[25];

% figgg=figure;
% bbb=get(figgg,'Position');
% h_factor=bbb(3)/bbb(4);
% new_width=8.7;
% helper_length=length(Aa);
% % load('Built_mat_discretePLoopaUP_1_time_0_85_amp.mat');
% % Median_mat_discrete(2)=Median_mat_discrete(2)/10;
% set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width length(Aa)/2*new_width]);

for kk=1:1:helper_length
% subplot(helper_length,1,kk);
non_eq_str = extract_non_eq_string(Rafiq_of_the_many{kk,:});

parts = strsplit(Rafiq_of_the_many{kk,:}, '\');

% Get the last part
last_part = parts{end};

TN=horzcat(num2str(Aa(kk)),'_T_',num2str(Aa2(kk)),'_N');
try
Adressme=horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\mimi_SS_god_',TN,'_',non_eq_str,'.mat');
load(Adressme);
Adressme1=horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\mimi_SS_god_discrete_',TN,'_',non_eq_str,'.mat');
load(Adressme1);
load(horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\mimi_SS_god_discrete_',TN,'_',non_eq_str,'.mat'));
catch
Adressme=horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\mimi_SS_god_',TN,non_eq_str,'_eq.mat');
load(Adressme);
Adressme1=horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\mimi_SS_god_discrete_',TN,non_eq_str,'_eq.mat');
load(Adressme1);
load(horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\mimi_SS_god_discrete_',TN,non_eq_str,'_eq.mat'));

end


% Adressme1=horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\mimi_SS_god_discrete_',TN,'_',non_eq_str,'.mat');
% load(Adressme1);

%try
%Adressme2=horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\Built_mat_discretePLoopaUP_Eq_',regexprep(TN,'N','P'),'.mat');
%load(Adressme2);
%catch
% Adressme2=horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\Built_mat_discretePLoopaUP',regexprep(TN,'N','P'),'.mat');
Adressme2=horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\Built_mat_discreteP',last_part,'.mat');
load(Adressme2);
%end


[X, Y] = extract_X_Y_N(Adressme);
X=num2str(X);
Y=num2str(Y);
%for eq refrence
% load(horzcat('energies_val',X,'_T_',Y,'_P.mat')); %give you energy_val
% energies_val=energies_true;

% Construct the file names dynamically using X and Y
built_mat_file = horzcat('Built_mat_discretePLoopaUP_N_Control_C_', X, '_T_', Y, '_P.mat');% gives Built_mat_discreteP
median_mat_file = horzcat('Median_mat_continuousLoopaUP_N_Control_C_', X, '_T_', Y, '_P.mat');% gives  Median_2
discrete_median_mat_file = horzcat('Median_mat_discreteLoopaUP_N_Control_C_', X, '_T_', Y, '_P.mat');% gives  Median_2
load(discrete_median_mat_file);

mimi_file = horzcat('mimi_SS_god_', X, '_T_', Y, '_N_eq'); %gives in eq mimi god
EQ_Built_mat_discreteP=load(built_mat_file);
EQ_Built_mat_discreteP=EQ_Built_mat_discreteP.Built_mat_discreteP;

load(median_mat_file);

mimi_file_discrete = horzcat('mimi_SS_god_discrete_', X, '_T_', Y, '_N_eq'); %gives in eq mimi god
mimi_file_discrete=load(mimi_file_discrete);
EQ_mimi_SS_god_discrete=mimi_file_discrete.mimi_SS_god_discrete;
% EQ_discrete_median_mat_file=load(discrete_median_mat_file);
% EQ_discrete_median_mat_file=discrete_median_mat_file.Median_mat_discrete;
ZZ_where_eq=squeeze(EQ_mimi_SS_god_discrete);

EQM=load(mimi_file);
EQ_mimi_SS_god=EQM.mimi_SS_god;
ZZ2=squeeze(EQ_mimi_SS_god);
if size(ZZ_where_eq,1)==length(energies_true)
energies_val=energies_true;
end
% load(horzcat('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\mimi_SS_god_discrete_',TN,'_',non_eq_str,'.mat'));
J=energies_true;
ZZ=squeeze(mimi_SS_god);
% ZZ_where_eq=squeeze(mimi_SS_god_discrete);
% ZZ(2,:)=ZZ(2,:)/10;
% maxZ_vec=zeros(1,length(J));

Mmm=zeros(1,size(ZZ,1));
Mmm_drive=zeros(1,size(ZZ,1));
ZZ_where_eq_yes_drive=squeeze(mimi_SS_god_discrete);
EQ_Built_mat_discreteP2=zeros(1,size(ZZ,1));
for yy=1:1:(length(J)) 
% maxZ=max(ZZ(yy,:));
% maxZ_vec=maxZ;
I=find(energies_val-J(yy)==0);
% maxZ=median(ZZ2(I,:));
Zaz=ZZ2(I,find(ZZ_where_eq(I,:)~=2*10^7));
Mmm(yy)=median(Zaz);

ZaZa=ZZ(yy,find(ZZ_where_eq_yes_drive(yy,:)~=2*10^7));
Mmm_drive(yy)=mean(ZaZa);
% ZZ(yy,find(ZZ_where_eq(yy,:)==2*10^7))=maxZ;
% ZZ(yy,:)=ZZ(yy,:)./maxZ;
EQ_Built_mat_discreteP2(yy)=EQ_Built_mat_discreteP(I);
end

% Mmm=median(ZZ,2);

end


% Mmm(2)=Mmm(2)/10;

% % figgg=figure;
% yyaxis left
% % bbb=get(figgg,'Position');
% % h_factor=bbb(3)/bbb(4);
% % new_width=8.7;
% % load('Built_mat_discretePLoopaUP_1_time_0_85_amp.mat');
% % Median_mat_discrete(2)=Median_mat_discrete(2)/10;
% % set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);
% %energies_true=[2 2.2 2.4 2.5 2.6 2.8 3 3.1 3.2 3.3 3.4 3.5];% 3.2 3.4
% % plot(flip(-energies_true,2),flip(Median_mat_discrete,1),'.b');
% % hold on;
% % plot(flip(-energies_true,2),flip(Mmm,1),'.r');
% % ylabel('Discrete Med(Tfas)')
% % xlabel('Js[k_{B}T]')
% % print (figgg,'-dpng','-r300');
% % 
% % aa=subplot(2,1,1);
% pop=Mmm;
% % figure; scatter (J,flip(ZZ),60,'filled');
% ylabel('{\boldmath${\hat{T}_{FAS}}$}','FontSize',6,'interpreter','latex')
% J=energies_true;
% % xlabel('$J_{s}$','FontSize',6,'interpreter','latex')
% hold on;
% plot(-flip(J,2),flip(pop',2),'Color',  [0.5 0.5 1],'LineWidth',2,'LineStyle','--');
% hold on;
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

% for yy=1:1:(length(J)) 
%     a=(ZZ(yy,:));
% [C,ia,ic] = unique(a);
% a_counts = accumarray(ic,1);
% try
% value_counts = [C, a_counts];
% catch
% value_counts = [C, a_counts'];
% 
% end
% % BubblePlot(ones(1,size(a_counts,1))*J(length(J)-yy+1), C,0.5*log(a_counts)+1,'b',8);
% dd=colormap(cool);%following gili request for xomplemntary
% bubbleplotUP(-ones(1,size(a_counts,1))*J(yy), C,0.5*log(a_counts)+1,dd(1,:),8,1,sum(a_counts),max(a_counts),size(ZZ,2));
% hold on;
% % BubblePlot(ones(1,size(1,1))*J(yy), C(end),0.5*log(a_counts(end))+1,'k',8);
% % hold on;
% end
% if kk==1
% cb = colorbar('northoutside');
% cbscalein = cb.Limits;
% cbscaleout = [0 100];
% cb.Ticks=0:0.2:1.1;
% cb.TickLabels = [{'0%','20%','40%','60%','80%','100%'}];
% end
% % title('{\boldmath${\hat{T}_{FAS}(Js),\, Eq}$}','FontSize',6,'interpreter','latex')
%  % ylim([0 105]);
% set(gca,'FontSize',6)
% hold on;
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
% yyaxis right
% % set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);
% %energies_true=[2 2.2 2.4 2.5 2.6 2.8 3 3.1 3.2 3.3 3.4 3.5];% 3.2 3.4
% plot(flip(-energies_true,2),100*flip(Built_mat_discreteP,1),'ok-','MarkerSize',4);
% % load('Built_mat_discretePLoopaUP_1_time_0_85_amp2.mat');
% % hold on;
% % plot(flip(-energies_true,2),100*flip(Built_mat_discreteP,1),'.r--','MarkerSize',6,'Marker','square');
% ylabel('{\boldmath${SA[\%]}$}','FontSize',6,'interpreter','latex');
% title(horzcat('{\boldmath${M_{T}=',num2str(Aa(kk)),',N=',num2str(Aa2(kk)),'}$}'),'FontSize',6,'interpreter','latex');
% if kk==helper_length
% xlabel('{\boldmath${J_{s}[k_{B}T]}$}','FontSize',6,'interpreter','latex');
% end
% ylim([-10 105]);
% set(gca,'YTick',0:20:100);

% ax=gca;
% ax.YAxis(1).Color = [0.5 0.5 1];
% ax.YAxis(2).Color = 'k';

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
% hold on;
% end
% cb.Position=[0.1311    0.9691    0.7744    0.0082];
% % set(gca,'FontSize',6);
% DOGG='Fig2EXTSI';
% cd('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC');
% print (DOGG,'-dpng','-r300');
% cd('C:\Users\admin\Pictures');
% print (DOGG,'-dpng','-r300');