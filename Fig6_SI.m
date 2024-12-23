current_dir=cd;
%assumed activation from the main directory, where this script is found in
addpath(horzcat(current_dir,'\Data\SI_Fig6'));


gg=figure; 
  bbb=get(gg,'Position');
  h_factor=bbb(3)/bbb(4);
  new_width=1*8.7;
  set(gg, 'Units', 'centimeters', 'Position',[2 2 new_width 1.5*new_width]);
ggx=gg.Number;
for kk=1:1:2
if kk==1
load('mimi_SS_god_2_T_25_N_eq.mat');
load('mimi_SS_god_discrete_2_T_25_N_eq.mat');
else
load('mimi_SS_god_2_T_25_N_non_eq_1_t_1_5_A.mat');
load('mimi_SS_god_discrete_2_T_25_N_non_eq_1_t_1_5_A.mat');
end
Time_Factor=4*10^-7;
energies =  [2 2.2 2.4 2.5 2.6 2.8 3 3.1 3.2 3.3 3.4 3.5 3.6];
% energies =  [ 3.4 3.5 3.6];
mimi_SS_god=squeeze(mimi_SS_god).*Time_Factor;
mimi_SS_god_discrete=squeeze(mimi_SS_god_discrete);
I=find(mimi_SS_god_discrete==2*10^7);
mimi_SS_god(I)=10^9;
CUTOFF=48;
for ii=1:1:3
if size(mimi_SS_god,1)~=3
mimi=mimi_SS_god(ii+10,:);
hh=2*ii-1;
Plot_Histograms_KLD_KMC_EFI(mimi(1:CUTOFF),0,energies(ii+10),1,Time_Factor,hh,ggx);
else
hh=2*ii;

mimi=mimi_SS_god(ii,:);
Plot_Histograms_KLD_KMC_EFI(mimi(1:CUTOFF),0,energies(ii),1,Time_Factor,hh,ggx);    
end
end
end

DOGG='FigS6';
cd(current_dir);
print (DOGG,'-dpng','-r300');
% cd('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC');
% print (DOGG,'-dpng','-r300');
