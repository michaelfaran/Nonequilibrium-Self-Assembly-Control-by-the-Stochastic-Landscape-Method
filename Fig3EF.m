current_dir=cd;
%assumed activation from the main directory, where this script is found in
addpath(horzcat(current_dir,'\Data\Fig3EF'));

time_Factor=4*10^-7; %michael add, the protein diffusion time scaling
% cd('C:\Users\admin\Documents\RMFO\LoopaUP_1_time_0_85_amp\3_4\0\07_28_2024_12_51_07_Js_3_4_Mu_0_Y_Beast_1_time_0_85_ampY_control_folder_num_0_job_id_9018968\2\9');
load('tot_energy_mu_0.0_energy_0.0_run_num_1_total_num_target_2.mat');
Total_energy=foo;
load('beast_cp_mu_0.0_energy_0.0_run_num_1_total_num_target_2.mat');
CP=foo;
cp=CP(:,1);
End_time=max(CP(:,4))*time_Factor;
load('mini_times_vec_mu_0.0_energy_0.0_run_num_1_total_num_target_2.mat');
time_o=foo;
%Fig3 for the paper
bumm=2*10^6;
fcc=2*10^7/bumm;
IL=1:1:2*10^7;
ILL=1:1*fcc:2*10^7;
% time_o=pivot4.foo; 
time= time_o(:,1)*time_Factor;
Aop=cumsum(time);
CCp=Aop(end);
load('mini_summed_distance_UP_vec_mu_0.0_energy_0.0_run_num_1_total_num_target_2.mat');
cd(current_dir);

Distance=foo;
qqqqq=Distance;
Distance=Distance(IL,:);
Aop2=Aop(IL);

% Original vector
x = Aop2;

% Sort the vector to make handling duplicates easier
x_sorted = sort(x);

% Define a very small number
epsilon = 1e-10;

% Loop through the sorted vector to add a small number to duplicates
for i = 2:length(x_sorted)
    if x_sorted(i) == x_sorted(i-1)
        x_sorted(i) = x_sorted(i) + epsilon;
        epsilon = epsilon * 10; % Increase the small number for next duplicates
    end
end

Aop2=x_sorted;

Total_energy2=interp1(Aop2, Total_energy,0:CCp/bumm:(CCp-CCp/bumm));
%         Distance(:,1)
SSSS=interp1(Aop2, single(Distance(:,1)),0:CCp/bumm:(CCp-CCp/bumm));
SSSS2=interp1(Aop2, single(Distance(:,2)),0:CCp/bumm:(CCp-CCp/bumm));
timi=0:CCp/bumm:CCp;
Distance=[int32(SSSS); int32(SSSS2)]';        
Total_energy2(1)=0;
Distance(1,1)=qqqqq(1,1);
Distance(1,2)=qqqqq(1,2);

Aop22=Aop2(ILL);
Time=0:CCp/bumm:(CCp-CCp/bumm); % the factor 4 is because of how the vector is itnerpolated
% Total_energy=o.data;
% load('Total_energy2.mat');
% cp=sort(o.trend.cp(1:o.trend.ncp_median));
% trend=o.trend.slp;

gg=figure; 
yyaxis left
  bbb=get(gg,'Position');
  h_factor=bbb(3)/bbb(4);
  new_width=0.8*8.7;
  set(gg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);

if sum(isnan(cp))~=0
   cp=cp(1:(find (isnan(cp),1)-1));     
end   
cp=[1 ;cp];
%taking the first one as well
%     cp=[1;cp];
if ~isempty(find (cp==0))
cp=cp(find (cp~=0)); 
end     
% noisy_signal = awgn(Total_energy, snr, 'measured');
        plot(Time,Total_energy2);
        ylimits=ylim;
        hold on;
        % L=(find(Distance(:,1)==0,1))*10^3;
        % cps=cp(cp<(min(find(Distance(:,2)==0,1),find(Distance(:,1)==0,1)))*fcc);
        cps=cp;
        ylim([-120 -60]);
        ylimits=ylim;
        %for i = 1 : length(cps)
            %plot( [cp(i),cp(i)]*10000, get(gca,'Ylim'),'color',[0.8500 0.3250 0.0980]);

            %plot( [cp(i),cp(i)]*10000, get(gca,'Ylim'),'color',[0.8500 0.3250 0.0980]);
xxxxx=Aop2(ceil((cps)));
s=scatter(xxxxx,(ylimits(1)+3).*ones(1,length(xxxxx)),[],[0.8500 0.3250 0.0980],"|",'filled','MarkerEdgeColor',[0.8500 0.3250 0.0980]);
            s.SizeData = 20;
         %    s=scatter(Aop2(ceil((cp(i)))),ylimits(1)+3,[],[0.8500 0.3250 0.0980],"|",'filled','MarkerEdgeColor',[0.8500 0.3250 0.0980]);
         % s.SizeData = 20;


            %xlim([0 2000]*10000);
            %ylabel('$E$\ $[K_{B}T]$','Interpreter','latex');
             % ylabel( '{\boldmath${E$\ $[K_{B}T]}$}','Interpreter','latex')
         ylabel( '{\boldmath${E\ [k_{B}T]}$}','Interpreter','latex')

            % xx0=xlabel('Time [Cs]');
            xx0= xlabel( '{\boldmath${Time\ [s]}$}','Interpreter','latex');
%              xx0.Position=xx0.Position-[0 10 0];
            % ylim([min(Total_energy)-10 0])  
            xlim([0 Time(end)])
        %end      
 xlim([ 22500*time_Factor 23440*time_Factor+300*time_Factor]);
 % xlim([Time(ceil((cps(end-25)/fcc))) Time(2317)*fcc]);
set(gca,'FontSize',6);
yyaxis right
Trend=CP(:,7);
time_for_trend2=CP(:,3)*time_Factor;
time_for_trend2=time_for_trend2(time_for_trend2>0);
Trend=Trend(time_for_trend2>0);
yaks=find(time_for_trend2<(23440+8200)*time_Factor);
time_for_trend=time_for_trend2(yaks);
Trend=Trend(yaks);

times=time_for_trend; trend_values=Trend/time_Factor;time_durations = diff(time_for_trend);  % Difference between consecutive time points
% Initialize time axis and square wave
times=[0; times];
total_time = times(end);  % The total time span from the first to the last time point
time_axis = times(1):time_Factor:total_time;  % A finer time axis for smooth plotting
square_wave = zeros(1, length(time_axis));  % Initialize the square wave

% Fill the square wave based on trend values and time intervals
current_time_index = 1;
trend_values(269)=-0.0003/time_Factor;

for i = 1:length(trend_values)-1
    % Find the time indices that correspond to the current time interval
    time_indices = time_axis >= times(i) & time_axis < times(i+1);
    
    % Assign the trend value to the corresponding time interval
    square_wave(time_indices) = trend_values(i);
end
% The last segment: fill in the last trend value until the last time
% square_wave(time_axis >= times(end-1)) = trend_values(end);


plot(time_axis, square_wave,'k'); 
xlim([ 22500*time_Factor 23440*time_Factor+300*time_Factor]);
ylim([-3 0.5]/time_Factor);
pp=ylim;
ax=gca;
xgca=gca;
ax.YAxis(2).Color = 'k';  
set(gca, 'Position',[0.1300 0.2214 0.7338 0.7036]);
 ylabel( '{\boldmath${\ Trend\ [{k_{B}T}\cdot{s^{-1}}]}$}','Interpreter','latex')
xgca.Position= [0.1700+0.08 0.2214 0.7338-0.18 0.7036];
xgca.Units = 'pixels';
xgca.Position= [29.6000+10   25.3540+10  145.7720 72.1600];
% [42.0960 24.2650 145.7720 72.1600]
DOGG='Fig3_maxime';


print (DOGG,'-dpng','-r300');

%Subfigure after control activated


gg=figure; 
yyaxis left
  bbb=get(gg,'Position');
  h_factor=bbb(3)/bbb(4);
  new_width=8.7;
  set(gg, 'Units', 'centimeters', 'Position',[2 2 0.5*new_width 0.5*new_width]);

if sum(isnan(cp))~=0
   cp=cp(1:(find (isnan(cp),1)-1));     
end   
cp=[1 ;cp];
%taking the first one as well
%     cp=[1;cp];
if ~isempty(find (cp==0))
cp=cp(find (cp~=0)); 
end     
% noisy_signal = awgn(Total_energy, snr, 'measured');
        plot(Time,Total_energy2);
        ylimits=ylim;
        hold on;
        % L=(find(Distance(:,1)==0,1))*10^3;
        % cps=cp(cp<(min(find(Distance(:,2)==0,1),find(Distance(:,1)==0,1)))*fcc);
        cps=cp;
        ylim([-150 -20]);
        ylimits=ylim;
        for i = 1 : length(cps)
            %plot( [cp(i),cp(i)]*10000, get(gca,'Ylim'),'color',[0.8500 0.3250 0.0980]);
         s=scatter(Aop2(ceil((cp(i)))),ylimits(1)+7,[],[0.8500 0.3250 0.0980],"|",'filled','MarkerEdgeColor',[0.8500 0.3250 0.0980]);
         s.SizeData = 20;
            %xlim([0 2000]*10000);
            %ylabel('$E$\ $[K_{B}T]$','Interpreter','latex');
             % ylabel( '{\boldmath${E$\ $[K_{B}T]}$}','Interpreter','latex')
         ylabel( '{\boldmath${E\ [k_{B}T]}$}','Interpreter','latex')

            % xx0=xlabel('Time [Cs]');
            xx0= xlabel( '{\boldmath${Time\ [s]}$}','Interpreter','latex');
%              xx0.Position=xx0.Position-[0 10 0];
            % ylim([min(Total_energy)-10 0])  
            xlim([0 Time(end)*time_Factor])
        end      
 % xlim([Time(ceil((cps(end-25)/fcc))) Time(2317)*fcc]);
set(gca,'FontSize',6);
yyaxis right
Trend=CP(:,7);

[I,M]=find(Total_energy2<-135,1);
last_time=Time(M);
time_for_trend2=CP(:,3)*time_Factor;
time_for_trend2=time_for_trend2(time_for_trend2>0);
Trend=Trend(time_for_trend2>0);
yaks=find(time_for_trend2<last_time+2300);
time_for_trend=time_for_trend2(yaks);
Trend=Trend(yaks);

times=time_for_trend; trend_values=Trend/time_Factor;time_durations = diff(time_for_trend);  % Difference between consecutive time points
% Initialize time axis and square wave
times=[0; times];
total_time = times(end);  % The total time span from the first to the last time point
time_axis = times(1):time_Factor:total_time;  % A finer time axis for smooth plotting
square_wave = zeros(1, length(time_axis));  % Initialize the square wave

% Fill the square wave based on trend values and time intervals
current_time_index = 1;
trend_values(269)=-0.0003/time_Factor;
for i = 1:length(trend_values)
    % Find the time indices that correspond to the current time interval
    time_indices = time_axis >= times(i) & time_axis < times(i+1);
    
    % Assign the trend value to the corresponding time interval
    square_wave(time_indices) = trend_values(i);
end
% The last segment: fill in the last trend value until the last time
% square_wave(time_axis >= times(end-1)) = trend_values(end);

yyaxis right
plot(time_axis, square_wave,'k'); 
xlim([ 22500*time_Factor last_time+300*time_Factor]);
ylim([-3 0.5]/time_Factor);
pp=ylim;
ax=gca;
ax.YAxis(2).Color = 'k';  
set(gca, 'Position',[0.1300 0.2214 0.7338 0.7036]);
 ylabel( '{\boldmath${\ Trend\ [{k_{B}T}\cdot{s^{-1}}]}$}','Interpreter','latex')
ax.Position= [0.1700+0.08 0.2214 0.7338-0.18 0.7036];
DOGG='Fig3_maxime_after_shock';
print (DOGG,'-dpng','-r300');



gg=figure; 
yyaxis left
  bbb=get(gg,'Position');
  h_factor=bbb(3)/bbb(4);
  new_width=8.7;
  set(gg, 'Units', 'centimeters', 'Position',[2 2 1*new_width 0.5*new_width]);

if sum(isnan(cp))~=0
   cp=cp(1:(find (isnan(cp),1)-1));     
end   
cp=[1 ;cp];
%taking the first one as well
%     cp=[1;cp];
if ~isempty(find (cp==0))
cp=cp(find (cp~=0)); 
end     
% noisy_signal = awgn(Total_energy, snr, 'measured');
        plot(Time,Total_energy2);
        ylimits=ylim;
        hold on;
        % L=(find(Distance(:,1)==0,1))*10^3;
        % cps=cp(cp<(min(find(Distance(:,2)==0,1),find(Distance(:,1)==0,1)))*fcc);
        cps=cp;
        ylim([-150 0]);
        ylimits=ylim;
        for i = 1 : length(cps)
            %plot( [cp(i),cp(i)]*10000, get(gca,'Ylim'),'color',[0.8500 0.3250 0.0980]);
         s=scatter(Aop2(ceil((cp(i)))),ylimits(1)+7,[],[0.8500 0.3250 0.0980],"|",'filled','MarkerEdgeColor',[0.8500 0.3250 0.0980]);
         s.SizeData = 20;
            %xlim([0 2000]*10000);
            %ylabel('$E$\ $[K_{B}T]$','Interpreter','latex');
             % ylabel( '{\boldmath${E$\ $[K_{B}T]}$}','Interpreter','latex')
         ylabel( '{\boldmath${E\ [k_{B}T]}$}','Interpreter','latex')

            % xx0=xlabel('Time [Cs]');
            xx0= xlabel( '{\boldmath${Time\ [s]}$}','Interpreter','latex');
%              xx0.Position=xx0.Position-[0 10 0];
            % ylim([min(Total_energy)-10 0])  
            xlim([0 Time(end)])
        end      
  % xlim([ 0 last_time]);
 % xlim([Time(ceil((cps(end-25)/fcc))) Time(2317)*fcc]);
set(gca,'FontSize',6);
yyaxis right
Trend=CP(:,7);
time_for_trend2=CP(:,3)*time_Factor;
time_for_trend2=time_for_trend2(time_for_trend2>0);
Trend=Trend(time_for_trend2>0);
yaks=find(time_for_trend2<last_time+5300*time_Factor);
time_for_trend=time_for_trend2(yaks);
Trend=Trend(yaks);

times=time_for_trend; trend_values=Trend/time_Factor;time_durations = diff(time_for_trend);  % Difference between consecutive time points
% Initialize time axis and square wave
times=[0; times];
total_time = times(end);  % The total time span from the first to the last time point
time_axis = times(1):time_Factor:total_time;  % A finer time axis for smooth plotting
square_wave = zeros(1, length(time_axis));  % Initialize the square wave

% Fill the square wave based on trend values and time intervals
current_time_index = 1;
trend_values(269)=-0.0003/time_Factor;

for i = 1:length(trend_values)-1
    % Find the time indices that correspond to the current time interval
    time_indices = time_axis >= times(i) & time_axis < times(i+1);
    
    % Assign the trend value to the corresponding time interval
    square_wave(time_indices) = trend_values(i);
end
% The last segment: fill in the last trend value until the last time
% square_wave(time_axis >= times(end-1)) = trend_values(end);
xlim([ 0 last_time+5300*time_Factor]);
plot(time_axis, square_wave,'k'); 
% xlim([ 22500 23440]*time_Factor);
%  ylim([-8 1.5]/time_Factor);
% pp=ylim;
ax=gca;
xgca=gca;
ax.YAxis(2).Color = 'k';  
set(gca, 'Position',[0.1300 0.2214 0.7338 0.7036]);
 ylabel( '{\boldmath${\ Trend\ [{k_{B}T}\cdot{s^{-1}}]}$}','Interpreter','latex')
% xgca.Position= [0.1700+0.08 0.2214 0.7338-0.18 0.7036];
DOGG='Fig3_up;_to_assembly';


print (DOGG,'-dpng','-r300');

