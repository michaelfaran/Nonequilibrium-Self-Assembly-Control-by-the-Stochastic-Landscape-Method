% clear;
% load('Fig4WS.mat');
current_cd='C:\Users\admin\Documents\Third_paper\Transfer\GitHub_Folders';
Bast_Path_Results=horzcat(current_cd,'\Data\Fig4');
cd(Bast_Path_Results);
Yo_dog='Fig4WS.mat';
load(Yo_dog);
cd(current_cd);
% Xn=2;
% Yn=36;
figgg=figure;
bbb=get(figgg,'Position');
h_factor=bbb(3)/bbb(4);
new_width=8.7;
% load('Built_mat_discretePLoopaUP_1_time_0_85_amp.mat');
% Median_mat_discrete(2)=Median_mat_discrete(2)/10;
set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 1*new_width]);
hold on;
subplot(2,1,1);
num_columns=size(Built_mat_combined,2);
% Define a colormap for different colors
colors = lines(num_columns);  % Get distinct colors for the plot

% Loop over each column and plot it with a different color
if search_me==1
[I,kp]=sort(t_values);
if length(kp)==15
    kp=kp(1:end-1);
end
Built_mat_combined=Built_mat_combined(kp,:);
if Xn==2 & Yn==25
Built_mat_combined=[EQ_Built_mat_discreteP; Built_mat_combined];
end
doubleMatrix = cell2mat(Mmm_drive_all');
doubleMatrix(isnan(doubleMatrix))=0;
doubleMatrix=doubleMatrix(kp,:);
doubleMatrix=[Mmm;doubleMatrix];
doubleMatrix(find(Built_mat_combined<0))=0;
t_values2=t_values(kp);

if Xn==2 & Yn==25
t_values=[1 t_values];
t_values2=t_values;
end
Built_mat_combined2=100*Built_mat_combined;
end
%yyaxis left   
for col = 1:num_columns
    % Plot the entire line (including all points) first for legend and connectivity
    h = plot(t_values2, Built_mat_combined2(:, col), 'Color', colors(col, :), 'LineWidth', 1,'LineStyle','-', 'Marker', 'none');
  hold on;  
    % Plot the first value with an empty circle marker (without affecting the legend)
    plot(t_values2(1), Built_mat_combined2(1, col), 'Color', colors(col, :), 'LineWidth', 1,'LineStyle','-', ...
        'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'none', 'HandleVisibility', 'off'); % Hide from legend
hold on;
    % Plot the rest of the values with filled circle markers (without affecting the legend)
    plot(t_values2(2:end), Built_mat_combined2(2:end, col), 'Color', colors(col, :), 'LineWidth', 1,'LineStyle','-', ...
        'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', colors(col, :), 'HandleVisibility', 'off'); % Hide from legend
hold on;
end
% Add legend, with the energies corresponding to each column
% hh=legend(arrayfun(@(x) sprintf('\\boldmath$Js=-%.1f$', x), energies, 'UniformOutput', false), ...
%        'FontSize', 6, 'Interpreter', 'latex', 'Location', 'best');
% hh.ItemTokenSize(1)=10;
% Label the axes
ylabel('{\boldmath${SA[\%]}$}','FontSize',6,'interpreter','latex');
% xlabel('{\boldmath${\rho}$}','FontSize',6,'interpreter','latex');
% title('Plot of Built_mat_discreteP columns vs. t values');
ylim([0 100]);
yticks((0:0.2:1)*100);
xlim([1 1.7]);
hh=legend(arrayfun(@(x) sprintf('\\boldmath$Js=-%.1f [k_{B}T]$', x), energies, 'UniformOutput', false), ...
       'FontSize', 6,'Location','northwest', 'Interpreter', 'latex');
hh.ItemTokenSize(1)=10;
% set(hh,"Position",[0.503433361319927,0.616895639641863,0.26661758702243,0.10304877917941]);
set(hh,"Location",'northwest');
if Xn==2 & Yn==25 
DOGG='Fig4';
set(hh,"Position",[0.503433361319927,0.616895639641863,0.26661758702243,0.10304877917941]);
elseif Xn==3 & Yn==25
set(hh,"Position",[0.403433361319927,0.616895639641863,0.26661758702243,0.10304877917941]); 
elseif Xn==4 & Yn==25
set(hh,"Position",[0.303433361319927,0.616895639641863,0.26661758702243,0.10304877917941]); 

else
    DOGG=horzcat('Fig4_',num2str(Yn),'P_',num2str(Xn),'_T');
end

set(gca, 'XTickLabel', []);
set(gca,'FontSize',6);
% print (DOGG,'-dpng','-r300');
% cd('C:\Users\admin\Pictures');
% print (DOGG,'-dpng','-r300');
% 
% % ax=gca;
% %  ax.YColor = 'k';
% % yyaxis right
% 
% figgg=figure;
% bbb=get(figgg,'Position');
% h_factor=bbb(3)/bbb(4);
% new_width=8.7;
text(-0.1, 1.1, 'A', 'Units', 'normalized', 'FontSize', 18, 'FontName', 'Calibri', 'HorizontalAlignment', 'center');

subplot(2,1,2);
% load('Built_mat_discretePLoopaUP_1_time_0_85_amp.mat');
% Median_mat_discrete(2)=Median_mat_discrete(2)/10;
% set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);
% hold on;


%plot(t_values,doubleMatrix)
doubleMatrix=doubleMatrix.*Time_Factor;
doubleMatrix= doubleMatrix./doubleMatrix(1,:);
for col = 1:num_columns
    % Plot the entire line (including all points) first for legend and connectivity
    h = plot(t_values, doubleMatrix(:, col), 'Color', colors(col, :), 'LineWidth', 1,  'LineStyle','-','Marker', 'none', 'HandleVisibility', 'on');
  hold on;  
    % Plot the first value with an empty circle marker (without affecting the legend)
    plot(t_values(1), doubleMatrix(1, col), 'Color', colors(col, :), 'LineWidth', 1,'LineStyle','-', ...
        'Marker', 's', 'MarkerSize', 4, 'MarkerFaceColor', 'none', 'HandleVisibility', 'off'); % Hide from legend
hold on;
    % Plot the rest of the values with filled circle markers (without affecting the legend)
    plot(t_values(2:end), doubleMatrix(2:end, col), 'Color', colors(col, :), 'LineWidth', 1,'LineStyle','-', ...
        'Marker', 's', 'MarkerSize', 4, 'MarkerFaceColor', colors(col, :), 'HandleVisibility', 'off'); % Hide from legend
hold on;
end
hold off;
% ax=gca;
%  ax.YColor = 'k';
ylabel('{\boldmath${\hat{T}_{FAS}}$}','FontSize',6,'interpreter','latex')
 % ylabel('{\boldmath${\hat{T}^_{FAS}[s]}$}','FontSize',6,'interpreter','latex');
ylim([0 1.5]);
yticks([0:0.5:1.5]);
xlabel('{\boldmath${\rho}$}','FontSize',6,'interpreter','latex');
xlim([1 1.7]);
text(-0.1, 1.1, 'B', 'Units', 'normalized', 'FontSize', 18, 'FontName', 'Calibri', 'HorizontalAlignment', 'center');

% Add legend, with the energies corresponding to each column
% legend(arrayfun(@(x) sprintf('\\boldmath$Js=-%.1f$', x), energies, 'UniformOutput', false), ...
%        'FontSize', 6, 'Interpreter', 'latex', 'Orientation', 'horizontal', 'Location', 'northoutside');
% hh=legend(arrayfun(@(x) sprintf('\\boldmath$Js=-%.1f$', x), energies, 'UniformOutput', false), ...
%        'FontSize', 6, 'Interpreter', 'latex', 'Location', 'northeast');
% hh.ItemTokenSize(1)=10;
% if Xn==3 | Xn==4
%     hh=legend(arrayfun(@(x) sprintf('\\boldmath$Js=-%.1f$', x), energies, 'UniformOutput', false), ...
%            'FontSize', 6, 'Interpreter', 'latex', 'Location', 'southwest');
%     hh.ItemTokenSize(1)=10;
% else
%     hh=legend(arrayfun(@(x) sprintf('\\boldmath$Js=-%.1f$', x), energies, 'UniformOutput', false), ...
%            'FontSize', 6, 'Interpreter', 'latex', 'Location', 'northeast');
%     hh.ItemTokenSize(1)=10;    
% end
set(gca,'FontSize',6);
if Xn==2 & Yn==25
DOGG='Fig4T';
else
    DOGG=horzcat('Fig4Tfas_',num2str(Yn),'P_',num2str(Xn),'_T');
end
print (DOGG,'-dpng','-r300');
% cd('C:\Users\admin\Pictures');
% print (DOGG,'-dpng','-r300');
