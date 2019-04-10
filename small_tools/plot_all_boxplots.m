function plot_all_boxplots()
% this function need read in all two instrument and EWS paired data table,
% then will try to make a boxplot that have all tables show up
size_fig = 1/2;
save_fig = 1;
save_records = 1;
plot_path = 'C:\Users\ZhaoX\Documents\paper\CF\fig\6-8\test\';
DMP_filter = false; % if use DMP filter, measurements with vortex abouve Eureka will be removed!
%SZA_filter = true; % if true, only use SZA <=85 (summer data type)
summer_spring_filter = 'false'; % if 'summer', only use SZA <=85 (summer data type), if 'spring', only use SAZ >85 (spring data type), or if 'false', all data will be used (no SZA filter)
use_SAOZ_V3_NDACC_SZArange = true; % if use SAOZ V3 processed with NDACC SZA range (86-91)
remove_SAOZ2017 = true;
remove_GBS2017_spring = true;

if save_records == true
    fid = fopen([plot_path 'relative_difference.txt'],'wt');
    fwrite(fid,'Relative Difference between X and Y /n');
    fwrite(fid,' Rel_diff = sum((X-Y)/((X+Y)/2))/N*100 (%) /n/n');
end
%addpath('E:\F\Work\MatlabCode');
addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab');
%load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews');
%load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews_2017'); % this version of data extended to 2017
%load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews_2017_high_quality_dmp.mat'); % this version of data extended to 2017
%load('C:\Users\ZhaoX\Documents\paper\CF\fig_old\gbs_saoz_brewer_merra2_ews_2017_high_quality_dmp.mat');
load('C:\Users\ZhaoX\Documents\paper\CF\fig\gbs_saoz_brewer_merra2_ews_2017_high_quality_dmp_nogbs2017spring.mat');
if use_SAOZ_V3_NDACC_SZArange == true
    load('C:\Projects\SAOZ\SAOZ_V3_NDACC_SZA_reformat_MERRA2_Brewer_EWS.mat');
end

SAOZ_V3_reformat_MERRA2_EWS(isnan(SAOZ_V3_reformat_MERRA2_EWS.mean_vcd),:) = [];
SAOZ_V3_reformat_Brewer_EWS(isnan(SAOZ_V3_reformat_Brewer_EWS.mean_vcd),:) = [];
SAOZ_V3_reformat_Brewer_EWS(isnan(SAOZ_V3_reformat_Brewer_EWS.SAOZ_V3_mean_vcd),:) = [];


if remove_SAOZ2017 == true
    TF_year = SAOZ_V3_reformat_Brewer_EWS.SAOZ_V3_UTC_str.Year == 2017;
    SAOZ_V3_reformat_Brewer_EWS(TF_year,:) = [];
    TF_year = SAOZ_V3_reformat_MERRA2_EWS.UTC_str.Year == 2017;
    SAOZ_V3_reformat_MERRA2_EWS(TF_year,:) = [];
end

if remove_GBS2017_spring == true
    GBS_MERRA2_EWS.Datetime = datetime(GBS_MERRA2_EWS.UTC_str);
    TF_2017spring = (GBS_MERRA2_EWS.Datetime < datetime('2017-04-30')) & (GBS_MERRA2_EWS.Datetime > datetime('2017-01-01'));
    GBS_MERRA2_EWS(TF_2017spring,:) = [];
    
    GBS_CF_MERRA2_EWS.Datetime = datetime(GBS_CF_MERRA2_EWS.UTC_str);
    TF_2017spring = (GBS_CF_MERRA2_EWS.Datetime < datetime('2017-04-30')) & (GBS_CF_MERRA2_EWS.Datetime > datetime('2017-01-01'));
    GBS_CF_MERRA2_EWS(TF_2017spring,:) = [];
    
    GBS_Brewer_EWS.Datetime = datetime(GBS_Brewer_EWS.UTC_str);
    TF_2017spring = (GBS_Brewer_EWS.Datetime < datetime('2017-04-30')) & (GBS_Brewer_EWS.Datetime > datetime('2017-01-01'));
    GBS_Brewer_EWS(TF_2017spring,:) = [];
    
    GBS_CF_Brewer_EWS.Datetime = datetime(GBS_CF_Brewer_EWS.UTC_str);
    TF_2017spring = (GBS_CF_Brewer_EWS.Datetime < datetime('2017-04-30')) & (GBS_CF_Brewer_EWS.Datetime > datetime('2017-01-01'));
    GBS_CF_Brewer_EWS(TF_2017spring,:) = [];
end
    

if DMP_filter
    TF_vortex = GBS_MERRA2_EWS.MERRA2_sPV_at_Theta490 >= 1.4e-4; % if we want remove measurements within vortex
    GBS_MERRA2_EWS(TF_vortex,:) = [];
    TF_vortex = SAOZ_MERRA2_EWS.MERRA2_sPV_at_Theta490 >= 1.4e-4; % if we want remove measurements within vortex
    SAOZ_MERRA2_EWS(TF_vortex,:) = [];

    TF_vortex = GBS_CF_MERRA2_EWS.MERRA2_sPV_at_Theta490 >= 1.4e-4; % if we want remove measurements within vortex
    GBS_CF_MERRA2_EWS(TF_vortex,:) = [];
    TF_vortex = SAOZ_CF_MERRA2_EWS.MERRA2_sPV_at_Theta490 >= 1.4e-4; % if we want remove measurements within vortex
    SAOZ_CF_MERRA2_EWS(TF_vortex,:) = [];
    
    TF_vortex = SAOZ_V3_reformat_MERRA2_EWS.MERRA2_sPV_at_Theta490 >= 1.4e-4; % if we want remove measurements within vortex
    SAOZ_V3_reformat_MERRA2_EWS(TF_vortex,:) = [];
   
end


%%
x = GBS_Brewer_EWS.mean_vcd./DU;
y = GBS_Brewer_EWS.mean_ColumnO3;
N=size(x);
%rel_diff = sum((x-y)./((x+y)/2))*100/N(1);
%rel_diff = sum((x-y)./(y))*100/N(1);
rel_diff = mean(x-y)/mean(y)*100;
disp(['GBS vs Brewer rel diff = ' num2str(rel_diff)]);

x = SAOZ_Brewer_EWS.mean_vcd./DU;
y = SAOZ_Brewer_EWS.mean_ColumnO3;
N=size(x);
%rel_diff = sum((x-y)./((x+y)/2))*100/N(1);
%rel_diff = sum((x-y)./(y))*100/N(1);
rel_diff = mean(x-y)/mean(y)*100;
disp(['SAOZ vs Brewer rel diff = ' num2str(rel_diff)]);


%% use MERRA-2 as benchmark
fig1 = figure;
hold all;
fig2 = figure;
hold all;

for i =1:5
    figure(fig1);% let's plot on fig1
    if i==1
        Table3 = GBS_MERRA2_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'MERRA2_Ozone';
        dataname = 'UT-GBS';
        %color =[0.5, 0.5, 0.5];
        color =[0.0, 0.45, 0.74];
    elseif i ==2
        Table3 = GBS_CF_MERRA2_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'MERRA2_Ozone';
        dataname = 'UT-GBS_CF';
        %color =[0.95, 0.5, 0.5];     
        color =[0.85, 0.33, 0.1];
    elseif i ==3
        Table3 = SAOZ_MERRA2_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'MERRA2_Ozone';
        dataname = 'SAOZ';
        %color =[0.5, 0.95, 0.5];
        color =[0.93, 0.69, 0.13];
    elseif i ==4
        Table3 = SAOZ_CF_MERRA2_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'MERRA2_Ozone';
        dataname = 'SAOZ_CF';
        %color =[0.5, 0.5, 0.95];
        color =[0.47, 0.67, 0.19];
    elseif i ==5
        Table3 = SAOZ_V3_reformat_MERRA2_EWS;
        tco_x = 'mean_vcd';
        tco_y = 'MERRA2_Ozone';
        dataname = 'SAOZ_V3';
        color =[0.95, 0.5, 0.95];
    end
    
   if i ~=5 % only use SZA filter for non SAOZ V3 data
%         if SZA_filter == true % the filter used to select only "summer" data
%             TF = Table3.sza_max >85;
%             Table3(TF,:) = [];
%         end
        if strcmp(summer_spring_filter,'summer')
            TF = Table3.sza_max >85;
            Table3(TF,:) = [];
        elseif strcmp(summer_spring_filter,'spring')
            TF = Table3.sza_max >85;
            Table3(~TF,:) = [];
        elseif strcmp(summer_spring_filter,'false')
        end
    end
        
        
    %TF = strcmp(Table3.weather_median_ampm,'Clear') | strcmp(Table3.weather_median_ampm,'Mainly Clear') | strcmp(Table3.weather_median_ampm,'Mostly Cloudy') | strcmp(Table3.weather_median_ampm,'Cloudy') | strcmp(Table3.weather_median_ampm,'Ice Crystals') | strcmp(Table3.weather_median_ampm,'Rain') | strcmp(Table3.weather_median_ampm,'Snow');
    TF = strcmp(Table3.weather_median_ampm,'Clear') | strcmp(Table3.weather_median_ampm,'Mainly Clear') | strcmp(Table3.weather_median_ampm,'Mostly Cloudy') | strcmp(Table3.weather_median_ampm,'Cloudy') | strcmp(Table3.weather_median_ampm,'Ice Crystals');
    Table3(~TF,:) = [];

    %x_location = 1:7;
    x_location = 1:5;
    eval(['Table3.delta_o3 = (Table3.' tco_x '- Table3.MERRA2_Ozone);']);
    %eval(['Table3.p_delta_o3 = 100.*(Table3.' tco_x '- Table3.MERRA2_Ozone)./Table3.' tco_y ';']);
    eval(['Table3.p_delta_o3 = 100.*(Table3.' tco_x '- Table3.MERRA2_Ozone)./((Table3.' tco_y '+ Table3.' tco_x ')/2);']);
    Table3 = sortrows(Table3,'weather_median_ampm','ascend');
    Table3(isnan(Table3.p_delta_o3),:)=[];
    y = mean(Table3.p_delta_o3);
    N_total_meas = height(Table3);
    y_err = std(Table3.p_delta_o3)/(N_total_meas.^0.5);
    disp(' ');
    disp('+++++++ All weather conditions ++++++');
    disp([dataname  ' vs.  MERRA-2 --> rel diff: ' num2str(y) ' +/- ' num2str(y_err)]);
    disp('------ by weather type -------');
    
    fprintf(fid,' \n');
    fprintf(fid,'+++++++ All weather conditions ++++++\n');
    fprintf(fid,[dataname  ' vs.  MERRA-2 --> rel diff: ' num2str(y) ' +/- ' num2str(y_err) '\n']);
    fprintf(fid,'------ by weather type ------- \n');
  
    
    mean_delta_o3_table = grpstats(Table3,'weather_median_ampm',{'mean','std','sem','numel'},'DataVars',{'MERRA2_Ozone','delta_o3','p_delta_o3'});
    mean_delta_o3_table = sortrows(mean_delta_o3_table,'weather_median_ampm');
    all_weathers = unique(Table3.weather_median_ampm);
    N_weathers = numel(all_weathers);
    x = 1:N_weathers;
    %y = 100.*mean_delta_o3_table.mean_delta_o3./mean_delta_o3_table.mean_MERRA2_Ozone;
    %y_err = 100.*mean_delta_o3_table.sem_delta_o3./mean_delta_o3_table.mean_MERRA2_Ozone;
    y = mean_delta_o3_table.mean_p_delta_o3;
    y_err = mean_delta_o3_table.std_p_delta_o3./(mean_delta_o3_table.GroupCount.^0.5);
    disp([dataname ' mean bias with MERRA-2 ----']);
    for i_weather = 1:N_weathers
        disp([all_weathers{i_weather} ' : ' num2str(y(i_weather)) ' +/- ' num2str(y_err(i_weather))]);
        fprintf(fid,[all_weathers{i_weather} ' : ' num2str(y(i_weather)) ' +/- ' num2str(y_err(i_weather)) '\n']);
    end
    bar(x_location+(i-1)/6,y,'FaceColor',color,'BarWidth',1/6,'EdgeColor',[1 1 1]);
    %errorbar(x_location+(i-1)/6,y,y_err,'.','Color',color);
    errorbar(x_location+(i-1)/6,y,y_err,'.','Color','k');
    
    boxplot(Table3.p_delta_o3, Table3.weather_median_ampm,'positions',x_location+(i-1)/6,'Widths',1/6,'Colors',color,'Symbol','','MedianStyle','target','Whisker',0);
    
    figure(fig2); % let's plot on fig2
    bar(x_location+(i-1)/6,mean_delta_o3_table.numel_delta_o3,'FaceColor',color,'BarWidth',1/6,'EdgeColor',[1 1 1]);
    
end

figure(fig1); % let's finish fig1 with labels and legend
ylim([-10 10]);
ylabel('TCO % difference (X - MERRA-2) [%]');
xlabel('EWS reported weather');
%legend('GBS','','GBS_C_F','','SAOZ','','SAOZ_C_F','','SAOZ_V_3');
grid on;
rotateXLabels(gca,45);
print_setting(size_fig,save_fig,['MERRA-2_vs_GBSandSAOZ_boxplots']);

figure(fig2); % let's finish fig2 with labels and legend
ylim([0 610]);
ylabel('Number of coincident measurements with MERRA-2');
xticklabels(unique(Table3.weather_median_ampm));
xlabel('EWS reported weather');
legend('UT-GBS','UT-GBS_C_F','SAOZ','SAOZ_C_F','SAOZ_V_3');
grid on;
rotateXLabels(gca,45);
print_setting(size_fig,save_fig,['MERRA-2_vs_GBSandSAOZ_coincidentnumber']);

%% use Brewer as benchmark
fig1 = figure;
hold all;
fig2 = figure;
hold all;

for i =1:5
    figure(fig1);% let's plot on fig1
    if i==1
        Table3 = GBS_Brewer_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'mean_ColumnO3';
        dataname = 'UT-GBS';
        %color =[0.5, 0.5, 0.5];
        color =[0.0, 0.45, 0.74];
    elseif i ==2
        Table3 = GBS_CF_Brewer_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'mean_ColumnO3';
        dataname = 'UT-GBS_CF';
        %color =[0.95, 0.5, 0.5];  
        color =[0.85, 0.33, 0.1];
    elseif i ==3
        Table3 = SAOZ_Brewer_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'mean_ColumnO3';
        dataname = 'SAOZ';
        %color =[0.5, 0.95, 0.5];
        color =[0.93, 0.69, 0.13];
    elseif i ==4
        Table3 = SAOZ_CF_Brewer_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'mean_ColumnO3';
        dataname = 'SAOZ_CF';
        %color =[0.5, 0.5, 0.95];
        color =[0.47, 0.67, 0.19];
    elseif i ==5
        Table3 = SAOZ_V3_reformat_Brewer_EWS;
        tco_x = 'SAOZ_V3_mean_vcd';
        tco_y = 'mean_ColumnO3';
        dataname = 'SAOZ_V3';
        color =[0.95, 0.5, 0.95];
    end
    
    if i ~= 5 % only use SZA filter for non SAOZ V3 data
%         if SZA_filter == true % the filter used to select only "summer" data
%             TF = Table3.sza_max >85;
%             Table3(TF,:) = [];
%         end
        if strcmp(summer_spring_filter,'summer')
            TF = Table3.sza_max >85;
            Table3(TF,:) = [];
        elseif strcmp(summer_spring_filter,'spring')
            TF = Table3.sza_max >85;
            Table3(~TF,:) = [];
        elseif strcmp(summer_spring_filter,'false')
        end
    end
    
    %TF = strcmp(Table3.weather_median_ampm,'Clear') | strcmp(Table3.weather_median_ampm,'Mainly Clear') | strcmp(Table3.weather_median_ampm,'Mostly Cloudy') | strcmp(Table3.weather_median_ampm,'Cloudy') | strcmp(Table3.weather_median_ampm,'Ice Crystals') | strcmp(Table3.weather_median_ampm,'Rain') | strcmp(Table3.weather_median_ampm,'Snow');
    TF = strcmp(Table3.weather_median_ampm,'Clear') | strcmp(Table3.weather_median_ampm,'Mainly Clear') | strcmp(Table3.weather_median_ampm,'Mostly Cloudy') | strcmp(Table3.weather_median_ampm,'Cloudy') | strcmp(Table3.weather_median_ampm,'Ice Crystals');
    Table3(~TF,:) = [];

    %x_location = 1:7;
    x_location = 1:5;
    eval(['Table3.delta_o3 = (Table3.' tco_x '- Table3.mean_ColumnO3);']);
    %eval(['Table3.p_delta_o3 = 100.*(Table3.' tco_x '- Table3.mean_ColumnO3)./Table3.' tco_y ';']);
    eval(['Table3.p_delta_o3 = 100.*(Table3.' tco_x '- Table3.mean_ColumnO3)./((Table3.' tco_y '+ Table3.' tco_x ')/2);']);
    Table3 = sortrows(Table3,'weather_median_ampm','ascend');
    Table3(isnan(Table3.p_delta_o3),:)=[];
    y = mean(Table3.p_delta_o3);
    N_total_meas = height(Table3);
    y_err = std(Table3.p_delta_o3)/(N_total_meas.^0.5);
    disp(' ');
    disp('+++++++ All weather conditions ++++++');
    disp([dataname  ' vs.  Brewer --> rel diff: ' num2str(y) ' +/- ' num2str(y_err)]);
    disp('------ by weather type -------');
    
    fprintf(fid,' \n');
    fprintf(fid,'+++++++ All weather conditions ++++++\n');
    fprintf(fid,[dataname  ' vs.  MERRA-2 --> rel diff: ' num2str(y) ' +/- ' num2str(y_err) '\n']);
    fprintf(fid,'------ by weather type ------- \n');
    
    mean_delta_o3_table = grpstats(Table3,'weather_median_ampm',{'mean','std','sem','numel'},'DataVars',{'mean_ColumnO3','delta_o3','p_delta_o3'});
    mean_delta_o3_table = sortrows(mean_delta_o3_table,'weather_median_ampm');
    all_weathers = unique(Table3.weather_median_ampm);
    N_weathers = numel(all_weathers);
    x = 1:N_weathers;
    %y = 100.*mean_delta_o3_table.mean_delta_o3./mean_delta_o3_table.mean_mean_ColumnO3;
    %y_err = 100.*mean_delta_o3_table.sem_delta_o3./mean_delta_o3_table.mean_mean_ColumnO3;
    y = mean_delta_o3_table.mean_p_delta_o3;
    y_err = mean_delta_o3_table.std_p_delta_o3./(mean_delta_o3_table.GroupCount.^0.5);
    disp([dataname ' mean bias with Brewer ----']);
    for i_weather = 1:N_weathers
        disp([all_weathers{i_weather} ' : ' num2str(y(i_weather)) ' +/- ' num2str(y_err(i_weather))]);
    end
    bar(x_location+(i-1)/6,y,'FaceColor',color,'BarWidth',1/6,'EdgeColor',[1 1 1]);
    %errorbar(x_location+(i-1)/6,y,y_err,'.','Color',color);
    errorbar(x_location+(i-1)/6,y,y_err,'.','Color','k');
    
    boxplot(Table3.p_delta_o3, Table3.weather_median_ampm,'positions',x_location+(i-1)/6,'Widths',1/6,'Colors',color,'Symbol','','MedianStyle','target','Whisker',0);
    
    figure(fig2); % let's plot on fig2
    bar(x_location+(i-1)/6,mean_delta_o3_table.numel_mean_ColumnO3,'FaceColor',color,'BarWidth',1/6,'EdgeColor',[1 1 1]);
    
end
figure(fig1); % let's finish fig1 with labels and legend
ylim([-10 10]);
ylabel('TCO % difference (X - Brewer) [%]');
xlabel('EWS reported weather');
%legend('GBS','','GBS_C_F','','SAOZ','','SAOZ_C_F','','SAOZ_V_3');
grid on;
rotateXLabels(gca,45);
print_setting(size_fig,save_fig,['Brewer_vs_GBSandSAOZ_boxplots']);

figure(fig2); % let's finish fig2 with labels and legend
ylim([0 610]);
ylabel('Number of coincident measurements with Brewer');
xticklabels(unique(Table3.weather_median_ampm));
xlabel('EWS reported weather');
legend('UT-GBS','UT-GBS_C_F','SAOZ','SAOZ_C_F','SAOZ_V_3');
grid on;
rotateXLabels(gca,45);
print_setting(size_fig,save_fig,['Brewer_vs_GBSandSAOZ_coincidentnumber']);