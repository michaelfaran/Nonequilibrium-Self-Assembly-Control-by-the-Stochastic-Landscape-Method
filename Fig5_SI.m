
%This code creates the matricies of the trapping regions for all drives and

%Assuming cd is in GitHub_Folders 
current_cd=cd;
cd_address=cd;
% cd_address=cd('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC');
Bast_Path_Results=horzcat(current_cd,'\Data\Beast_Results');
% Bast_Path_Results='C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\Beast_Map_Results';
% RESULTS_PATH='C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC\Beast_Results_Many\LoopaUP_N_Control_C_2_T_36_P';
RESULTS_PATH='C:\Users\admin\Documents\Third_paper\Transfer\GitHub_Folders\Data\Fig3D\No_Control_C_2_T_25_P';
time_Factor=4*10^-7;
energies=[2 2.2 2.4 2.5 2.6 2.8 3 3.1 3.2 3.3 3.4 3.5 3.6];
energies=[3.6];
energies_val=energies;
% energies='2_2';
% drives=string(9);
drives=0:0.2:2.8;
drives=0;
drives_val=drives;

Matricia1=zeros(length(energies),length(drives));
Matricia2=zeros(length(energies),length(drives));
Matricia1_Factor=zeros(length(energies),length(drives));

for kk=1:1:length(energies)
for gg=1:1:length(drives)
mu=drives(gg);
energy=energies(kk);

cd(RESULTS_PATH);
filenameX=horzcat('C_Js_',regexprep(num2str(energy, '%5.1f'),'\.','_'),'_mu_',regexprep(num2str(mu, '%5.1f'),'\.','_'));
C2=load(filenameX);
cd(cd_address);


C_reduced=C2.D_Reduced;


%Here there can be both 5000 and 0 traj, lets see if they are
%distiniugshable
ssz=size(C_reduced);
CCC=C_reduced;

% bottom_E=-220;
% Upper_E=-168;
 Lower_t=-10*10^-3;
 Upper_t=10*10^-3;


%create a heat map of 2D
% figure;
x=CCC(:,5);
y=CCC(:,7);
C=CCC(:,8);

[i1,i2]=find(y>0);
[S,i4]=sort(y(i1));
miki=C(i1);
Sa=miki(i4);
window_size=0.002;
numba= round((max(y(i1))-mod(max(y(i1)),window_size))/window_size);

vecta=zeros(1,numba);
Xvecta=zeros(1,numba);
for ii=1:1:(numba+1)
imi=(ii-1)*window_size;
boost=window_size;
vecta(ii)=median(Sa(find(imi<S & S<imi+boost)));
Xvecta(ii)=imi;
end

lord_vecta= interpolate_nans(vecta);

[i1,i2]=find(y<0);
[S,i4]=sort(y(i1));
miki=C(i1);
Sa=miki(i4);
window_size=0.002;
numba_neg= round((min(y(i1))-mod(min(y(i1)),window_size))/window_size);

vecta_neg=zeros(1,numba_neg);
Xvecta_neg=zeros(1,numba_neg);
for ii=1:1:(abs(numba_neg)+1)
imi=(ii-1)*window_size;
boost=window_size;
vecta_neg(ii)=median(Sa(find(-imi-boost<S & S<-imi)));
Xvecta_neg(ii)=-imi;
end
lord_vecta_neg= interpolate_nans(vecta_neg);

VECTA=[lord_vecta_neg lord_vecta];
X_VECTA=[Xvecta_neg Xvecta];
gog=flip(0:0.01:1,2);
range_prob=0.2;
for lr=1:1:length(gog)
top_line=Xvecta(find(lord_vecta<gog(lr)*max(VECTA),1));
bottom_line=Xvecta_neg(find(lord_vecta_neg<gog(lr)*max(VECTA),1));
Lower_t=bottom_line;
Upper_t=top_line;
Precenta=100*(length(C(find(y<Upper_t & y>Lower_t))))/(length(C(find(y>Upper_t | y<Lower_t)))+length(C(find(y<Upper_t & y>Lower_t))));
Facta=median(C(find(y<Upper_t & y>Lower_t)))/median(C(find(y<Upper_t | y>Lower_t)));

if Precenta>20
    sprintf('The percentage of data insiade the traping region is:')
    display(Precenta)
    sprintf('Its dwelling time average is this factor smaller:')
    display(Facta);

break
end
end
figg=figure;
bbb=get(figg,'Position');
h_factor=bbb(3)/bbb(4);
new_width=8.7;
set(figg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.8*new_width]);
scatter(y/time_Factor,C.*time_Factor,1,'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
% xlabel('trend value');
% ylabel('Dwelling Time');
% xx0=xlabel('Time [Cs]');
% xx0= xlabel( '{\boldmath${Time\ [Cs]}$}','Interpreter','latex');
xx0=xlabel( '{\boldmath${Trend\ [{k_{B}T}\cdot{s^{-1}}]}$}','Interpreter','latex');
ylabel( '{\boldmath${Dwelling \ Time\ [s]}$}','Interpreter','latex');

hold on;
nameX='1D_map_of_dwell_time_vs_trend';
nameXTotal=horzcat(nameX,'_', regexprep(num2str(energy, '%5.1f'),'\.','_'),'_',regexprep(num2str(mu, '%5.1f'),'\.','_'));
xline(Lower_t/time_Factor,'k--','LineWidth',1);
xline(Upper_t/time_Factor,'k--','LineWidth',1);
% cd(Bast_Path_Results);
xlim ([-0.5 0.5]/time_Factor);
ylim([-2 1200]*time_Factor);
set(gca,'FontSize',6);
set(gca, 'box', 'on') 
cd(current_cd);
print (nameXTotal,'-dpng','-r300');

cd(cd_address);
Matricia1(kk,gg)=Lower_t;
Matricia2(kk,gg)=Upper_t;
Matricia1_Factor(kk,gg)=Facta;
close(figg);
end
end

Matricia_bottom_trend=Matricia1;
Matricia_top_trend=Matricia2;
% cd('C:\Users\admin\Documents\GitHub\GradProject')
% save('Matricia_bottom_trend','Matricia_bottom_trend');
% save('Matricia_top_trend','Matricia_top_trend');
% save('Matricia1_Factor','Matricia1_Factor');
% save('drives_val','drives_val');
% save('energies_val','energies_val');
% 
% cd(cd_address);

% sz=40;
% scatter(x,y,sz,C,'filled');
% hold on;
% colormap jet;
% colorbar;
% %caxis([100000 300000]);
% caxis([0 2*median(C)]);
% xlabel('M*');
% a=ylabel('t*');
% set(a,'rotation',0,'VerticalAlignment','middle');
% title('Macrostate Dwell Time');
% set(gca,'FontSize',16);
% yline(Lower_t,'k--','LineWidth',3);
% yline(Upper_t,'k--','LineWidth',3);
%  % ylim([30*Lower_t 30*Upper_t]);
% cd(Bast_Path_Results);
% %for 3_5, choosing the trend as 0.15 is better than 0.0002 
% print ('Time_vs_E_and_Trend_for_slides','-dpng','-r600');
% set(gcf,'color','w');
% %create a heat map of 2D
% figure;
% x=CCC(:,4);
% y=CCC(:,3);
% C=CCC(:,1);
% sz=40;
% scatter(x,y,sz,C,'filled');
% hold on;
% colormap jet;
% colorbar;
% caxis([min(C) max(C)]);
% xlabel('Time');
% ylabel('Trend');
% title('E');
% print ('E_vs_Time_and_Trend','-dpng','-r600');
% %create a heat map of 2D
% figure;
% x=CCC(:,1);
% y=CCC(:,4);
% C=CCC(:,3);
% sz=40;
% scatter(x,y,sz,C,'filled');
% hold on;
% colormap jet;
% colorbar;
% caxis([-0.4 max(C)]);
% xlabel('E');
% ylabel('Time');
% title('Trend');
% print ('Trend_vs_Time_and_E','-dpng','-r600');
% 
% x=CCC(:,1);
% y=CCC(:,3);
% C=CCC(:,4);
% gg=0.45;
% A=find_symmetric_ylines_simple(y,gg);
% y1=A(1);
% y2=A(2);
% figure;
% plot(y,C,'b.');
% yx=(abs(y1)+abs(y2))/2;
% hold on; xline(-yx);
% hold on;xline(yx);
% 
% 
% % II=size(I);
% % min_length=10000;
% % 
% % for ii=1:1:(II(1)-1)
% %  min_length=min(length(C_reduced(I(ii)+1:I(ii+1),1)),min_length);
% % end
% % 
% % indes_min=find(I(2:end)-I(1:(end-1))==min_length );
% % Ix1=I(indes_min(1));
% % Ix2=Ix1+min_length;
% %     x=C_reduced(Ix1+1:Ix2,1);%mean_vec
% %     y=C_reduced(Ix1+1:Ix2,2);%std_vec
% %     z=C_reduced(Ix1+1:Ix2,3);%trend
% %     c=C_reduced(Ix1+1:Ix2,6);
% % TrajectoryPlotVecs(x',y',z',c,cmin,cmax);
% % 
% % 
% 
% % xlabel('Mean','FontSize',28);
% % ylabel('Variance','FontSize',28);
% % zlabel('Trend','FontSize',28);
% % title('NEW COORD,Trajectories of Simulation CG States with Time Remaining to FIRST Self Assemble Colorbar')
% % % savefig('Stochastic 3');
% % savefig(figure(gcf),fullfile(Main_dir_hazirim,folderName,horzcat('Stochastic_6_', regexprep(num2str(mu, '%5.1f'),'\.','_'),'.fig')),'compact');


