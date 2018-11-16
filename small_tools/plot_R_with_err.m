function plot_R_with_err(output)
%% this is a function to plot scatter fitting parameters with error bars
% this function need to read in the output from "output = plot_all_scatterplots()"

save_fig = 0;
fitting_parameter = 'R';
load('scatter_status.mat');
output = output_sixweather;
column_nms = output.Properties.VariableNames;
row_nms = output.Properties.RowNames;
figure;hold all;
for i =1:numel(column_nms)
    column_nm = cell2mat(column_nms(i));
    for j = 1:numel(row_nms)
        row_nm = cell2mat(row_nms(j));
        eval([fitting_parameter '(i,j) = output.' column_nm '({''' row_nm '''},:).' fitting_parameter ';']);
        eval(['RU(i,j) = output.' column_nm '({''' row_nm '''},:).RU;']);
        eval(['RL(i,j) = output.' column_nm '({''' row_nm '''},:).RL;']);
    end
    start_position_factor = 1/(numel(column_nms)+2);
    x = (start_position_factor*i):numel(row_nms);
    eval(['y =' fitting_parameter]);
    %plot(x,y(i,:),'.-');
    neg = (R - RL);
    pos = (RU - R);
    errorbar(x,y(i,:),neg(i,:),pos(i,:));
end
xticks(x);
%xticklabels({'Clear','Cloudy','Mainly Clear','Mostly Cloudy','Ice Crystals','Rain','Snow'});
set(gca,'XTickLabel',row_nms)
legend_str = strrep(column_nms,'_','-');
legend(legend_str);
ylabel(fitting_parameter);

print_setting(1/2,save_fig,[fitting_parameter]);
