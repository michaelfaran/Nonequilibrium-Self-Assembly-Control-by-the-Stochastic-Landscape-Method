%This is the code that generates the figure 5 of the article
Acid=[1.1,1.4,1.7];
Hallmark_cd=cd;
 boxy=['A';'B';'C'];
  boxy2=['D';'E';'F'];
    CR = 3600;
figgg=figure;
bbb=get(figgg,'Position');
h_factor=bbb(3)/bbb(4);
new_width=17.8;
% load('Built_mat_discretePLoopaUP_1_time_0_85_amp.mat');
% Median_mat_discrete(2)=Median_mat_discrete(2)/10;
set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 1*new_width]);
t=tiledlayout(3,2);

for gg=1:1:length(Acid)
       Acid_Drive=Acid(gg);
       cd(horzcat(Hallmark_cd,'\Data\Fig5')) 
       load(horzcat('Homeboy_Js_3.4_Acid_',num2str(Acid_Drive),'_mu_0_0.mat'))   ;
       cd(Hallmark_cd);
        % Convert each cell array to double matrices
        double_before = cell2mat(d_mean_before_array);
        double_after = cell2mat(d_mean_after_array);
        double_after_5 = cell2mat(d_mean_after_5_array);
        double_after_7 = cell2mat(d_mean_after_7_array);
        double_after_9 = cell2mat(d_mean_after_9_array);

        double_during = cell2mat(d_mean_during_array);
        
        % Reshape each matrix to [numCells x 2] format
        double_before = reshape(double_before, [], 2);
        double_after = reshape(double_after, [], 2);
        double_after_5 = reshape(double_after_5, [], 2);
        double_after_7 = reshape(double_after_7, [], 2);
        double_after_9 = reshape(double_after_9, [], 2);
        double_during = reshape(double_during, [], 2);
        
        Duval=[min(double_after');min(double_after_5');min(double_after_7');min(double_after_9')];
                customColormap = [
            0.7, 0.9, 1.0;  % Light blue
            0.4, 0.7, 0.9;  % Medium-light blue
            0.2, 0.5, 0.8;  % Medium blue
            0.0, 0.3, 0.6;  % Dark blue
        ];
        
    A = min(double_before');
    C = min(double_during');


nexttile;
x = A;
III = find(x < CR);
x = x(III);
y = C(III);

% Define number of bins for averaging
numBins = 5;  % Adjust as necessary

% Group data into bins and calculate bin centers
[N, binEdges, binIndices] = histcounts(x, numBins);
binCenters = linspace(10, 90, numBins);  % Set bin centers manually for positioning
if isempty(find(N < 10))
    disp('N check');
end

% Calculate data points for each bin and store them in a cell array
binnedData = arrayfun(@(i) y(binIndices == i), 1:numBins, 'UniformOutput', false);

% Plot y = x line for reference
xRange = linspace(0, 100, 100);

pp=4;

% Plot individual boxplots at specified positions with wider boxes
for i = 1:numBins
    % Use 'boxplot' with specific x position for each bin center
    h=boxplot(binnedData{i}, 'Positions', binCenters(i), 'Widths', 15, 'Colors', customColormap(pp,:), 'Symbol', 'o');
outliers = findobj(h, 'Tag', 'Outliers');

% Change outlier marker properties to have a pale blue edge
for i = 1:length(outliers)
    set(outliers(i), 'Marker', 'o', 'MarkerSize', 6)
end
hold on;
end

p2=plot(xRange, xRange, '--k', 'LineWidth', 2);
    hold on;
% % Scatter plot of bin centers vs. bin averages
% yAverages = cellfun(@mean, binnedData);  % Calculate means for each bin
% scatter(binCenters, yAverages, 50, 'MarkerFaceColor', customColormap(pp,:), 'MarkerEdgeColor', customColormap(pp,:));

% Set x-axis ticks to range from 0 to 100 with ticks every 10 units
xlim([0 100]);
set(gca, 'XTick', 0:10:100, 'XTickLabel', 0:10:100);  % Manually set x-tick labels from 0 to 100

% Customize plot appearance
ylim([-5 120]);
if  gg==3  
xlabel('{ \boldmath$\langle d \rangle_{before}$}', 'Interpreter', 'latex');
end
ylabel('{ \boldmath$\langle d \rangle_{during}$}', 'Interpreter', 'latex');
% text(0.1, 0.9, boxy(gg), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold'); % Add label
text(-0.15, 1.05, boxy(gg), 'Units', 'normalized', 'FontSize', 18, 'FontName', 'Calibri', 'HorizontalAlignment', 'center');
title(['{\boldmath$\rho = ' num2str(Acid_Drive) '$}'], 'Interpreter', 'latex');

if gg~=3
set(gca, 'XTickLabel', []);
end

set(gca,'FontSize',6);

nexttile;
for pp = 4:1:size(Duval,1)

    B = Duval(pp,:);
    CR = 3600;
    
    x = A;
    III = find(x < CR);
    x = x(III);
    y = B(III);
    
    % Define number of bins for averaging
    numBins = 5;  % Adjust as necessary
    
    % Group data into bins and calculate bin centers
    [N1, binEdges, binIndices] = histcounts(x, numBins);
    if isempty(find(N1 < 10))
        disp('N1 check')
    end
    binCenters = linspace(10, 90, numBins);  % Manually set bin centers within the 0-100 range
    
    % Calculate data points for each bin and store them in a cell array
    binnedData = arrayfun(@(i) y(binIndices == i), 1:numBins, 'UniformOutput', false);

    % Plot y = x line for reference
    xRange = linspace(0, 100, 100);

    
    % Plot individual boxplots at specified positions with wider boxes
    for i = 1:numBins
        % Use 'boxplot' with specific x position for each bin center
        h=boxplot(binnedData{i}, 'Positions', binCenters(i), 'Widths', 15, 'Colors', customColormap(pp,:), 'Symbol', 'o');
   % Find the outliers
outliers = findobj(h, 'Tag', 'Outliers');

% Change outlier marker properties to have a pale blue edge
for i = 1:length(outliers)
    set(outliers(i), 'Marker', 'o', 'MarkerSize', 6)
end    
hold on;
    end  

    plot(xRange, xRange, '--k', 'LineWidth', 2);
    hold on;
    % Scatter plot of bin centers vs. bin averages
    % yAverages = cellfun(@mean, binnedData);  % Calculate means for each bin
    % scatter(binCenters, yAverages, 50, 'MarkerFaceColor', customColormap(pp,:), 'MarkerEdgeColor', customColormap(pp,:));

    % Set x-axis ticks to range from 0 to 100 with ticks every 10 units
    xlim([0 100]);
    set(gca, 'XTick', 0:10:100, 'XTickLabel', 0:10:100);  % Manually set x-tick labels from 0 to 100
ylim([-5 120]); 
    % Customize plot appearance
if  gg==3   
xlabel('{ \boldmath$\langle d \rangle_{before}$}', 'Interpreter', 'latex');
end
ylabel('{ \boldmath$\langle d \rangle_{after}$}', 'Interpreter', 'latex');
    % title('d after vs. d before Box Plot Average');
    % legend('y = x', 'Bin Averages', 'Location', 'best');
    hold on;

end
% text(0.1, 0.9, boxy2(gg), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold'); % Add label
text(-0.15, 1.05, boxy2(gg), 'Units', 'normalized', 'FontSize', 18, 'FontName', 'Calibri', 'HorizontalAlignment', 'center');
title(['{\boldmath$\rho = ' num2str(Acid_Drive) '$}'], 'Interpreter', 'latex');

if gg~=3
set(gca, 'XTickLabel', []);
end
% filenameXx = horzcat('Correlation_between_before_after', regexprep(num2str(E, '%5.1f'), '\.', '_'), '_Acid_', regexprep(num2str(Acid_Drive, '%5.1f'), '\.', '_'), '_mu_', regexprep(num2str(mu(j), '%5.1f'), '\.', '_'));
% print(filenameXx, '-dpng', '-r300');

% title('d during vs. d before Box Plot');
% hold off;

% Save the figure
% filenameXx = horzcat('Correlation_between_before_during', regexprep(num2str(E, '%5.1f'), '\.', '_'), '_Acid_', regexprep(num2str(Acid_Drive, '%5.1f'), '\.', '_'), '_mu_', regexprep(num2str(mu(j), '%5.1f'), '\.', '_'));
% print(filenameXx, '-dpng', '-r300');
set(gca,'FontSize',6);
end
hold on;
dummyPlot = plot(nan, nan, '-', 'LineWidth', 1.5, 'Color', customColormap(pp,:)); % Mimics boxplot color
% legend('y = x', 'Bin Averages', 'Location', 'southeast');
legend([p2, dummyPlot], {'$\ \ \ y = x$', 'Boxplots'}, 'Location', 'southeast', 'Interpreter', 'latex');

% Save the figure
cd(Hallmark_cd);
filenameXx = horzcat('Fig5');
print(filenameXx, '-dpng', '-r300');
% cd('C:\Users\admin\Pictures');
% print(filenameXx, '-dpng', '-r300');
