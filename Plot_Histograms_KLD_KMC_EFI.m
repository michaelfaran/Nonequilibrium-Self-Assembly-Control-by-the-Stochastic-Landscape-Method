%Plot all the histogram etc
function Plot_Histograms_KLD_KMC(mimi,mu,energy,use_time_flag,Time_Factor,hh,gg)


pop=median(mimi);
%assuming that for a guassian it is abour 34% of the data far away to the
%std in each side
mimi_sorted=sort(mimi);
pop_minus=mimi_sorted(floor(0.16*length(mimi)));
pop_plus=mimi_sorted(floor(0.84*length(mimi)));
% aa=figure;
% close(aa)
%  upp=histogram(mimi,length(mimi));
% upp=histogram(mimi);
%  if max(upp.BinEdges)<5*10^7
%     edges=[upp.BinEdges 5*10^7];
%     edges=linspace(0,edges(end),length(edges));
%     upp=histogram(mimi,edges);
%  else   
%  edges=[upp.BinEdges(1:(end-1)) 5*10^7];
%  edges=linspace(0,edges(end),length(edges));
%  upp=histogram(mimi,edges);
%  hold on;
%   bar(upp.BinEdges(2:end)-(upp.BinEdges(2)-upp.BinEdges(1))/2,[zeros(length(edges)-2,1);upp.Values(end)],1,'FaceColor','g'	)
%      
%  end

%Use unequal binning method, reference from 
%lognormal values binning 
figure;
mimi_normal=log(mimi(mimi<10^9));
mean_n=mean(mimi_normal);
upp=histogram(mimi_normal);
binning=upp.BinEdges;
edgy=exp(binning);
edgy(1)=0;
if use_time_flag==1
edgy(end)=1*max(mimi(mimi<10^9));
else
edgy(end)=2*10^7;
end

end_of_days=edgy(end);
a11=[edgy(edgy<pop)];
a1= [a11 edgy(length(a11))];
close(gcf);
figure;
upp2=histogram(mimi(mimi<edgy(end)));
binning2=upp2.BinEdges;
edgy2=binning2;
rrr=sort(a1);
a2=edgy2(edgy2>rrr(end));
edgy=[rrr a2(2:end)];

if edgy(end)<end_of_days
    edgy=[edgy end_of_days-1*Time_Factor end_of_days+end_of_days/20];
else
    edgy=edgy(1:find(edgy==max(edgy),1));
    edgy=[edgy(1:(end-1)) end_of_days-1*Time_Factor  end_of_days+end_of_days/20];
end

if use_time_flag==1
mimi(mimi==10^9)=end_of_days+1*end_of_days/20;
end
 close(gcf);
% gg=figure; 
%   bbb=get(gg,'Position');
%   h_factor=bbb(3)/bbb(4);
%   new_width=1*8.7;
%   set(gg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);
figure(gg); 
subplot(3,2,hh);
upp=histogram(mimi(mimi< max(mimi)),edgy);

 % fullscreen();
% figure; scatter(mu,mimi(:),60,'filled');
ylabel('{\boldmath${Count}$}','Interpreter','latex','FontSize',6);
% xlabel('T_{FAS}','FontSize',6);
if hh==5 || hh==6 
xx0= xlabel( '{\boldmath${T_{FAS}\ [s]}$}','Interpreter','latex','FontSize',6);
title( '{\boldmath${J_{s}=-3.6\ [k_{B}T]}$}','Interpreter','latex','FontSize',6);
elseif hh==3 || hh==4 
title( '{\boldmath${J_{s}=-3.5\ [k_{B}T]}$}','Interpreter','latex','FontSize',6);
elseif hh==1 || hh==2
title( '{\boldmath${J_{s}=-3.4\ [k_{B}T]}$}','Interpreter','latex','FontSize',6);

end
%title(horzcat('Tfas Distribution of mu ' , num2str(mu),' Strong Energy ',num2str(energy), ' Number of Targets ',str2),'FontSize',36);
set(gca,'FontSize',6)
mimi_sorted=sort(mimi);
pop=median(mimi);
pop_plus=mimi_sorted(floor(0.84*length(mimi)));
normal_diturb_std_low=log(pop/pop_minus);
normal_diturb_std_high=log(pop_plus/pop);
hold on;
% scatter(pop,1,90,'d','k','filled');
[MM,iiii]=min(abs(upp.BinEdges-pop));
% try
% x= upp.BinEdges(iiii)+upp.BinWidth/2;
%  [M,II]=min(abs(upp.BinEdges-pop));
%  y=upp.Values(II);
% % bar(x,y,upp.BinWidth,'r');

linear_scale_std_minus=exp(mean_n)-exp(mean_n-normal_diturb_std_low);
% title(horzcat('Med_{N}=',num2str(log(pop),'%0.2e'),', \mu_{N}=',num2str(mean_n,'%0.2e'),' , \sigma_{N+}=',num2str(normal_diturb_std_high,'%0.2e'),' , \sigma_{N-}= ',num2str(normal_diturb_std_low,'%0.2e'),' , \sigma_{N-,l}= ',num2str(linear_scale_std_minus,'%0.2e')));
% catch
x= upp.BinEdges(end)-0.5*(upp.BinEdges(end)-upp.BinEdges(end-1));

if use_time_flag==1
y=length(find(mimi==max(mimi)));
x= upp.BinEdges(end)+0.5*(upp.BinEdges(end)-upp.BinEdges(end-1));
else
y=length(find(mimi==2*10^7));
end
%  [M,II]=min(abs(upp.BinEdges-pop));
%  y=upp.Values(II-1);
hold on;
 bar(x,y,upp.BinEdges(end)-upp.BinEdges(end-1),'g');
if pop==pop_plus
c1=xline(x,'--r','LineWidth',1);
c21=xline(x,'--m','LineWidth',1);
c21=xline(pop_minus,'--m','LineWidth',1); 
else
c1=xline(pop,'--r','LineWidth',1);
c21=xline(pop_plus,'--m','LineWidth',1);
c21=xline(pop_minus,'--m','LineWidth',1); 
end

he=x+upp.BinEdges(end)-upp.BinEdges(end-1);
he=0.5;
% Get current xticks

xlim([0 he]);
xticks_current = get(gca, 'XTick');

% Check if he is within the current xticks
if ~ismember(he, xticks_current)
    % Calculate the interval between xticks (assumes a constant interval)
    xtick_interval = xticks_current(2) - xticks_current(1);
    
    % Generate new xticks until surpassing 'he'
    new_xticks = xticks_current;
    while max(new_xticks) < he
        new_xticks = [new_xticks, max(new_xticks) + xtick_interval]; %#ok<AGROW> 
    end
    
    % Set the new xticks
    set(gca, 'XTick', new_xticks);
end
% Get the current xticks
xticks_current = get(gca, 'XTick');

% Convert xticks to a cell array of character strings
xticks_as_cell = cellstr(num2str(xticks_current', '%.1f'));  % Format as needed (e.g., 1 decimal place)
 set(gca, 'XTickLabel', xticks_as_cell);
 xlim([0 xticks_current(end)]);
% c1=xline(pop,'--r','LineWidth',6);
% end
%legend('Histogram Data','Median')
if use_time_flag==1
saveas(gcf,horzcat('Continuous_Tfas_KMC_Histogram_of_mu_' , num2str(mu,'%.1f'),'_Strong_Energy_',num2str(energy,'%.1f'),'.png'));
saveas(gcf,horzcat('Continuous_Tfas_KMC_Histogram_of_mu_' , num2str(mu,'%.1f'),'_Strong_Energy_',num2str(energy,'%.1f'),'.fig'));
else
saveas(gcf,horzcat('Discrete_Tfas_KMC_Histogram_of_mu_' , num2str(mu,'%.1f'),'_Strong_Energy_',num2str(energy,'%.1f'),'.png'));
saveas(gcf,horzcat('Discrete_Tfas_KMC_Histogram_of_mu_' , num2str(mu,'%.1f'),'_Strong_Energy_',num2str(energy,'%.1f'),'.fig'));
end


% close(gcf);
% 
% gh=figure;
% if use_time_flag==1
% mimi_2=mimi(mimi<max(mimi));
% else
% mimi_2=mimi(mimi<2*10^7);
% end
%cancel the idea of taking out the late samples, can we think of a rule of
%thumb when this is relveant? should we conider for example showing results
%only when the total counts inside the extraordinary bin is less then half o
% mimi_2=mimi;
% [f,xi] = ksdensity(mimi,'support','positive');
% plot(xi,f,'LineWidth',1.3);
% try
% pd = fitdist(mimi_2','lognormal');
% hold on;
%  y = pdf(pd,xi);
%  plot(xi,y,'LineStyle','--','LineWidth',1.3)
% % figure; scatter(mu,mimi(:),60,'filled');
% % ylabel('P','FontSize',36);
% end
%  fullscreen();
% try
% kde= horzcat('KDE (',' \mu = ' ,num2str(mean(pd),2),' , ','\sigma = ', num2str(std(pd),2),' ) ');
% end
% lognormal= horzcat('Lognormal Fit (',' \mu = ' ,num2str(mean(f.*xi*(xi(2)-xi(1)))*length(xi),2),' , ','\sigma = ', num2str(std(f.*xi*(xi(2)-xi(1)))*length(xi),2),' ) ');
% try
% legend(kde,lognormal)
% end
% xlabel('T_{FAS}','FontSize',36);
% % xlim([0 5*10^7])
% % title(horzcat('Tfas Distribution of mu ' , num2str(mu),' Strong Energy ',num2str(energy), ' Number of Targets ',str2),'FontSize',36);
% set(gca,'FontSize',36);
% if use_time_flag==1
% saveas(gcf,horzcat('Continuous_Tfas_KMC_Distribution_of_mu_' , num2str(mu,'%.1f'),'_Strong_Energy_',num2str(energy,'%.1f'),'.png'));
% a=regexprep(horzcat('Results_of_Continuous_Tfas_KMC_Distribution_of_mu_' , num2str(mu,'%.1f'),'_Strong_Energy_',num2str(energy,'%.1f')),'\.','_'); 
% A=convertStringsToChars(horzcat(a,'.mat'));
% close(gh)
% save(A,'mimi');
% else
% saveas(gcf,horzcat('Discrete_Tfas_KMC_Distribution_of_mu_' , num2str(mu,'%.1f'),'_Strong_Energy_',num2str(energy,'%.1f'),'.png'));
% a=regexprep(horzcat('Results_of_Discrete_Tfas_KMC_Distribution_of_mu_' , num2str(mu,'%.1f'),'_Strong_Energy_',num2str(energy,'%.1f')),'\.','_'); 
% A=convertStringsToChars(horzcat(a,'.mat'));
% close(gh)
% save(A,'mimi');
end



 

