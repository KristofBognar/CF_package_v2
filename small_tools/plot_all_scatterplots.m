function output = plot_all_scatterplots()
% this function is used to plot all scatterplots between differen pairs of
% measurements/model. It gives multi-scatter plots in different weather
% conditions
save_fig = 1;
DU = 2.6870e+16;
%plot_path = 'E:\H\work\Eureka\GBS\CI\archive\sactter_plots\temp\';
plot_path = 'C:\Users\ZhaoX\Documents\paper\CF\fig\8\correct_R_err\'
%load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews');
%load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews_2017_high_quality_dmp.mat');% this dataset extend to 2017, and single measurement has been filtered
load('C:\Users\ZhaoX\Documents\paper\CF\fig\gbs_saoz_brewer_merra2_ews_2017_high_quality_dmp.mat');
%addpath('E:\F\Work\MatlabCode');
addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab\');

cd(plot_path);
% make groups of weather types we are interested
weather_types = {'All','Clear','Cloudy','Mainly Clear','Mostly Cloudy','Ice Crystals','Rain','Snow'};
% make groups of instrument/model pairs we are interested
instrument_pairs = {'GBS_Brewer_EWS','GBS_CF_Brewer_EWS','SAOZ_Brewer_EWS','SAOZ_CF_Brewer_EWS', ...
    'SAOZ_V3_reformat_Brewer_EWS','GBS_MERRA2_EWS','GBS_CF_MERRA2_EWS','SAOZ_MERRA2_EWS',...
    'SAOZ_CF_MERRA2_EWS','SAOZ_V3_reformat_MERRA2_EWS'};
output = table;



for i = 1:numel(instrument_pairs)

    for j = 1:numel(weather_types)
        h_all(j) = subplot(4,2,j); % ax for subplot
    end
    
    for j = 1:numel(weather_types)
        
        eval(['data = ' cell2mat(instrument_pairs(i)) ';']); % assign the data we will use
        
        if ~strcmp('All',weather_types(j)) % if not intend to plot all weather condition together
            TF = strcmp(data.weather_median_ampm,weather_types(j));% let's filter data by weather type
            data(~TF,:) = [];
        end
        
        if contains( cell2mat(instrument_pairs(i)) , 'Brewer') % if we use Brewer as benchmark
            y = data.mean_ColumnO3; % this is Brewer ozone column
            y_label = 'Brewer';
        elseif contains( cell2mat(instrument_pairs(i)) , 'MERRA2') % if we use Brewer as benchmark
            y = data.MERRA2_Ozone; % this is MERRA-2 ozone column
            y_label = 'MERRA-2';
        end
        
        if contains( cell2mat(instrument_pairs(i)) , 'SAOZ_V3_reformat_Brewer') % if we want compare with SAOZ V3
            x = data.SAOZ_V3_mean_vcd; % this is SAOZ V3 ozone column when paired with Brewer
            x_label = 'SAOZ_V_3';
        elseif contains( cell2mat(instrument_pairs(i)) , 'SAOZ_V3_reformat_MERRA2') % if we want compare with SAOZ V3
            x = data.mean_vcd; % this is SAOZ V3 ozone column hen paired with MERRA-2
            x_label = 'SAOZ_V_3';
        else % if we want compare with our CF package processed data
            x = data.mean_vcd./DU; % this is GBS,GBS_CF,SAOZ,and SAOZ_CF ozone column
            if contains( cell2mat(instrument_pairs(i)) , 'GBS') % if we want compare with GBS data
                if contains( cell2mat(instrument_pairs(i)) , 'CF')
                    x_label = 'GBS_C_F';
                else
                    x_label = 'GBS';
                end
            elseif contains( cell2mat(instrument_pairs(i)) , 'SAOZ') % if we want compare with GBS data
                if contains( cell2mat(instrument_pairs(i)) , 'CF')
                    x_label = 'SAOZ_C_F';
                else
                    x_label = 'SAOZ';
                end
            end
        end
        
        % extra filter to make sure no NaN will be used!
        TF = isnan(x) | isnan(y);
        x(TF,:) = [];
        y(TF,:) = [];
        
        % make some scatter plots, and give outputs of R, Number of
        % measurements, intercepts, and slop
        [R,RL,RU,N,k_1,intercept_1,k_2] = generic_scatter(x,y,x_label,y_label,weather_types(j),save_fig);
        h1_scatter = gca;
        copyobj(allchild(h1_scatter),h_all(j));
        xlabel(h_all(j),[x_label '_.']);
        ylabel(h_all(j),[y_label]);
        title(h_all(j),weather_types(j));
        xlim(h_all(j),[200 600]);
        ylim(h_all(j),[200 600]);
        close(gcf);
        % prepare a structure to be saved in output table
        st.R = R; % correlation coef
        st.RL = RL; % Lower bound for correlation coefficient
        st.RU = RU; % Upper bound for correlation coefficient
        st.N = N; % number of measurements
        st.k_1 = k_1; % slop of "y = a*x + b" fitting
        st.intercept_1 = intercept_1; % intercept of "y = a*x + b" fitting
        st.k_2 = k_2; % slop of "y = a*x" fitting
        
        % save above fitting results into a single output table
        weather_type = cell2mat(weather_types(j));% weather type name
        %weather_type= weather_type(find(~isspace(weather_type))); % remove
        %space 
        instrument_pair = cell2mat(instrument_pairs(i));% instrument pair name
        %eval(['output.' weather_type '({''' instrument_pair '''},1)=
        %st;']); % make table use weather type as column and instrument
        %pair as row
        eval(['output.' instrument_pair '({''' weather_type  '''},1)= st;']);% make table use weather type as row and instrument
        %pair as column
    end
    print_setting(1,save_fig,[y_label '_vs_' x_label '_all_merged' ]);
    close all;
end
plot_scatter_status(output,'R',save_fig);
plot_scatter_status(output,'N',save_fig);
plot_scatter_status(output,'k_1',save_fig);
plot_scatter_status(output,'intercept_1',save_fig);
plot_scatter_status(output,'k_2',save_fig);
%% this is a generic scatter plots will be used for differen inputs
function [R,RL,RU,N,k_1,intercept_1,k_2] = generic_scatter(x,y,x_label,y_label,weather_type,save_fig)
save_fig = 0;
%[R,N,k_1,intercept_1,k_2] = linear_fits(x,y);
[intercept,slop,slop_nlm,mdl_lm,mdl_nlm] = line_fits(x,y);
[R,P,RL,RU] = corrcoef(x,y);
R = R(1,2);
RL = RL(1,2);
RU = RU(1,2);
N = numel(x);
k_1 = slop;
intercept_1 = intercept;
k_2 = slop_nlm;
xlim([200 600]);
ylim([200 600]);
xlabel([x_label ' TCO [DU]']);
ylabel([y_label ' TCO [DU]']);
title(weather_type);
grid on;
print_setting(1/4,save_fig,[y_label '_vs_' x_label '_' cell2mat(weather_type)]);


%% this is a function to plot scatter fitting parameters
function plot_scatter_status(output,fitting_parameter,save_fig)
column_nms = output.Properties.VariableNames;
row_nms = output.Properties.RowNames;
figure;hold all;
for i =1:numel(column_nms)
    column_nm = cell2mat(column_nms(i));
    for j = 1:numel(row_nms)
        row_nm = cell2mat(row_nms(j));
        eval([fitting_parameter '(i,j) = output.' column_nm '({''' row_nm '''},:).' fitting_parameter ';']);
    end
    x = 1:numel(row_nms);
    eval(['y =' fitting_parameter]);
    plot(x,y(i,:),'.-');
end
xticks(x);
%xticklabels({'Clear','Cloudy','Mainly Clear','Mostly Cloudy','Ice Crystals','Rain','Snow'});
set(gca,'XTickLabel',row_nms)
legend_str = strrep(column_nms,'_','-');
legend(legend_str);
ylabel(fitting_parameter);

print_setting(1/2,save_fig,[fitting_parameter]);



