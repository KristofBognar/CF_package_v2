function plot_all_boxplots()
% this function need read in all two instrument and EWS paired data table,
% then will try to make a boxplot that have all tables show up
size_fig = 1/2;
save_fig = 1;

load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews');
SAOZ_V3_reformat_MERRA2_EWS(isnan(SAOZ_V3_reformat_MERRA2_EWS.mean_vcd),:) = [];
SAOZ_V3_reformat_Brewer_EWS(isnan(SAOZ_V3_reformat_Brewer_EWS.mean_vcd),:) = [];

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
        color =[0.5, 0.5, 0.5];
    elseif i ==2
        Table3 = GBS_CF_MERRA2_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'MERRA2_Ozone';
        color =[0.95, 0.5, 0.5];        
    elseif i ==3
        Table3 = SAOZ_MERRA2_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'MERRA2_Ozone';
        color =[0.5, 0.95, 0.5];
    elseif i ==4
        Table3 = SAOZ_CF_MERRA2_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'MERRA2_Ozone';
        color =[0.5, 0.5, 0.95];
    elseif i ==5
        Table3 = SAOZ_V3_reformat_MERRA2_EWS;
        tco_x = 'mean_vcd';
        tco_y = 'MERRA2_Ozone';
        color =[0.95, 0.5, 0.95];
    end
    TF = strcmp(Table3.weather_median_ampm,'Clear') | strcmp(Table3.weather_median_ampm,'Mainly Clear') | strcmp(Table3.weather_median_ampm,'Mostly Cloudy') | strcmp(Table3.weather_median_ampm,'Cloudy') | strcmp(Table3.weather_median_ampm,'Ice Crystals') | strcmp(Table3.weather_median_ampm,'Rain') | strcmp(Table3.weather_median_ampm,'Snow');
    Table3(~TF,:) = [];

    x_location = 1:7;
    eval(['Table3.delta_o3 = (Table3.' tco_x '- Table3.MERRA2_Ozone);']);
    eval(['Table3.p_delta_o3 = 100.*(Table3.' tco_x '- Table3.MERRA2_Ozone)./Table3.' tco_y ';']);
    Table3 = sortrows(Table3,'weather_median_ampm','ascend');
    
    
    mean_delta_o3_table = grpstats(Table3,'weather_median_ampm',{'mean','std','sem','numel'},'DataVars',{'MERRA2_Ozone','delta_o3'});
    mean_delta_o3_table = sortrows(mean_delta_o3_table,'weather_median_ampm');
    N_weathers = numel(unique(Table3.weather_median_ampm));
    x = 1:N_weathers;
    y = 100.*mean_delta_o3_table.mean_delta_o3./mean_delta_o3_table.mean_MERRA2_Ozone;
    y_err = 100.*mean_delta_o3_table.sem_delta_o3./mean_delta_o3_table.mean_MERRA2_Ozone;
    bar(x_location+(i-1)/10,y,'FaceColor',color,'BarWidth',0.08,'EdgeColor',[1 1 1]);
    %errorbar(x_location+(i-1)/10,y,y_err,'.','Color',color);
    errorbar(x_location+(i-1)/10,y,y_err,'.','Color','k');
    
    boxplot(Table3.p_delta_o3, Table3.weather_median_ampm,'positions',x_location+(i-1)/10,'Widths',0.09,'Colors',color,'Symbol','','MedianStyle','target','Whisker',0);
    
    figure(fig2); % let's plot on fig2
    bar(x_location+(i-1)/10,mean_delta_o3_table.numel_MERRA2_Ozone,'FaceColor',color,'BarWidth',0.08,'EdgeColor',[1 1 1]);
    
end

figure(fig1); % let's finish fig1 with labels and legend
ylim([-10 8]);
ylabel('TCO % difference (X - MERRA-2) [%]');
xlabel('EWS reported weather');
%legend('GBS','','GBS_C_F','','SAOZ','','SAOZ_C_F','','SAOZ_V_3');
grid on;
rotateXLabels(gca,45);
print_setting(size_fig,save_fig,['MERRA-2_vs_GBSandSAOZ_boxplots']);

figure(fig2); % let's finish fig2 with labels and legend
ylim([0 500]);
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
        color =[0.5, 0.5, 0.5];
    elseif i ==2
        Table3 = GBS_CF_Brewer_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'mean_ColumnO3';
        color =[0.95, 0.5, 0.5];        
    elseif i ==3
        Table3 = SAOZ_Brewer_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'mean_ColumnO3';
        color =[0.5, 0.95, 0.5];
    elseif i ==4
        Table3 = SAOZ_CF_Brewer_EWS;
        tco_x = 'mean_vcd./DU';
        tco_y = 'mean_ColumnO3';
        color =[0.5, 0.5, 0.95];
    elseif i ==5
        Table3 = SAOZ_V3_reformat_Brewer_EWS;
        tco_x = 'SAOZ_V3_mean_vcd';
        tco_y = 'mean_ColumnO3';
        color =[0.95, 0.5, 0.95];
    end
    TF = strcmp(Table3.weather_median_ampm,'Clear') | strcmp(Table3.weather_median_ampm,'Mainly Clear') | strcmp(Table3.weather_median_ampm,'Mostly Cloudy') | strcmp(Table3.weather_median_ampm,'Cloudy') | strcmp(Table3.weather_median_ampm,'Ice Crystals') | strcmp(Table3.weather_median_ampm,'Rain') | strcmp(Table3.weather_median_ampm,'Snow');
    Table3(~TF,:) = [];

    x_location = 1:7;
    eval(['Table3.delta_o3 = (Table3.' tco_x '- Table3.mean_ColumnO3);']);
    eval(['Table3.p_delta_o3 = 100.*(Table3.' tco_x '- Table3.mean_ColumnO3)./Table3.' tco_y ';']);
    Table3 = sortrows(Table3,'weather_median_ampm','ascend');
    
    
    mean_delta_o3_table = grpstats(Table3,'weather_median_ampm',{'mean','std','sem','numel'},'DataVars',{'mean_ColumnO3','delta_o3'});
    mean_delta_o3_table = sortrows(mean_delta_o3_table,'weather_median_ampm');
    N_weathers = numel(unique(Table3.weather_median_ampm));
    x = 1:N_weathers;
    y = 100.*mean_delta_o3_table.mean_delta_o3./mean_delta_o3_table.mean_mean_ColumnO3;
    y_err = 100.*mean_delta_o3_table.sem_delta_o3./mean_delta_o3_table.mean_mean_ColumnO3;
    bar(x_location+(i-1)/10,y,'FaceColor',color,'BarWidth',0.08,'EdgeColor',[1 1 1]);
    %errorbar(x_location+(i-1)/10,y,y_err,'.','Color',color);
    errorbar(x_location+(i-1)/10,y,y_err,'.','Color','k');
    
    boxplot(Table3.p_delta_o3, Table3.weather_median_ampm,'positions',x_location+(i-1)/10,'Widths',0.09,'Colors',color,'Symbol','','MedianStyle','target','Whisker',0);
    
    figure(fig2); % let's plot on fig2
    bar(x_location+(i-1)/10,mean_delta_o3_table.numel_mean_ColumnO3,'FaceColor',color,'BarWidth',0.08,'EdgeColor',[1 1 1]);
    
end
figure(fig1); % let's finish fig1 with labels and legend
ylim([-10 8]);
ylabel('TCO % difference (X - Brewer) [%]');
xlabel('EWS reported weather');
%legend('GBS','','GBS_C_F','','SAOZ','','SAOZ_C_F','','SAOZ_V_3');
grid on;
rotateXLabels(gca,45);
print_setting(size_fig,save_fig,['Brewer_vs_GBSandSAOZ_boxplots']);

figure(fig2); % let's finish fig2 with labels and legend
ylim([0 500]);
ylabel('Number of coincident measurements with Brewer');
xticklabels(unique(Table3.weather_median_ampm));
xlabel('EWS reported weather');
legend('UT-GBS','UT-GBS_C_F','SAOZ','SAOZ_C_F','SAOZ_V_3');
grid on;
rotateXLabels(gca,45);
print_setting(size_fig,save_fig,['Brewer_vs_GBSandSAOZ_coincidentnumber']);