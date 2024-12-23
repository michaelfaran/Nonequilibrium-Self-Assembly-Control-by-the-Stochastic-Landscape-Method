%Fig3 for the paper
current_cd=cd;
close all;
cd(horzcat(current_cd,'\Data\Fig3BC'));
load('TREND_o.mat');
time_Factor=4*10^-7;
Time=o.time/4; % the factor 4 is because of how the vector is itnerpolated
Total_energy=o.data;
load('Total_energy2.mat');
cp=sort(o.trend.cp(1:o.trend.ncp_median));
trend=o.trend.slp;
cd(current_cd);

gg=figure; 
  bbb=get(gg,'Position');
  h_factor=bbb(3)/bbb(4);
  new_width=8.7;
  set(gg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);

yyaxis("left");
if sum(isnan(cp))~=0
   cp=cp(1:(find (isnan(cp),1)-1));     
end   
cp=[1 ;cp];
%taking the first one as well
%     cp=[1;cp];
if ~isempty(find (cp==0))
cp=cp(find (cp~=0)); 
end     

         plot(Time.*1000*time_Factor,Total_energy);
         ylimits=ylim;
        hold on;
        for i = 1 : length(cp)
            %plot( [cp(i),cp(i)]*10000, get(gca,'Ylim'),'color',[0.8500 0.3250 0.0980]);
         s=scatter(Time(cp(i))*1000*time_Factor,ylimits(1)+5,[],[0.8500 0.3250 0.0980],"|",'filled','MarkerEdgeColor',[0.8500 0.3250 0.0980]);
         s.SizeData = 20;
            %xlim([0 2000]*10000);
            %ylabel('$E$\ $[K_{B}T]$','Interpreter','latex');
             % ylabel( '{\boldmath${E$\ $[K_{B}T]}$}','Interpreter','latex')
         ylabel( '{\boldmath${E\ [k_{B}T]}$}','Interpreter','latex')

            % xx0=xlabel('Time [Cs]');
            xx0= xlabel( '{\boldmath${Time\ [s]}$}','Interpreter','latex');
%              xx0.Position=xx0.Position-[0 10 0];
            ylim([min(Total_energy)-10 0])
            xlim([0 Time(end)*1000*time_Factor])
        end 
yyaxis("right");

plot(Time.*10000*time_Factor,o.trend.slp/time_Factor,'k'); 
pp=ylim;
ax=gca;
ax.YAxis(2).Color = 'k';  
set(gca, 'Position',[0.1300 0.2214 0.7338 0.7036]);
 ylabel( '{\boldmath${Trend\ [{k_{B}T}\cdot{s^{-1}}]}$}','Interpreter','latex')
DOGG='Fig3_maxime';
set(gca,'FontSize',6);
print (DOGG,'-dpng','-r300');

%1.5015*10^7, 1.7818*10^7

gg=figure; 
  bbb=get(gg,'Position');
  h_factor=bbb(3)/bbb(4);
  new_width=0.67*8.7;
  set(gg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);

yyaxis("left");
if sum(isnan(cp))~=0
   cp=cp(1:(find (isnan(cp),1)-1));     
end   
cp=[1 ;cp];
%taking the first one as well
%     cp=[1;cp];
if ~isempty(find (cp==0))
cp=cp(find (cp~=0)); 
end     

         plot(Time.*1000.*time_Factor,Total_energy);
         ylimits=ylim;
        hold on;
        for i = 1 : length(cp)
            %plot( [cp(i),cp(i)]*10000, get(gca,'Ylim'),'color',[0.8500 0.3250 0.0980]);
         s=scatter(Time(cp(i))*1000*time_Factor,-146.5,[],[0.8500 0.3250 0.0980],"|",'filled','MarkerEdgeColor',[0.8500 0.3250 0.0980]);
         s.SizeData = 20;
            %xlim([0 2000]*10000);
            %ylabel('$E$\ $[K_{B}T]$','Interpreter','latex');
             % ylabel( '{\boldmath${E$\ $[K_{B}T]}$}','Interpreter','latex')
         % ylabel( '{\boldmath${E\ [K_{B}T]}$}','Interpreter','latex')

            % xx0=xlabel('Time [Cs]');
            % xx0= xlabel( '{\boldmath${Time\ [Cs]}$}','Interpreter','latex');
%              xx0.Position=xx0.Position-[0 10 0];
            ylim([min(Total_energy)-10 -100])
            xlim([0 Time(end)*1000*time_Factor])
        end 
                 ylabel( '{\boldmath${E\ [k_{B}T]}$}','Interpreter','latex')

yyaxis("right");
 ylabel( '{\boldmath${Trend\ [{k_{B}T}\cdot{s^{-1}}]}$}','Interpreter','latex')
            xx0= xlabel( '{\boldmath${Time\ [s]}$}','Interpreter','latex');

plot(Time.*1000*time_Factor,o.trend.slp/time_Factor,'k'); 
xlim([1.5015*10^6, 1.7818*10^6]*time_Factor);
ylim([-15 5]./time_Factor);
ax=gca;
ax.YAxis(2).Color = 'k';
set(gca,'FontSize',6);
ax.Position=[    0.1368+0.05    0.2115    0.7626-0.1    0.6560];
DOGG='Fig3_minime';
 print (DOGG,'-dpng','-r300');