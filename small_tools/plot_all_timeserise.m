function plot_all_timeserise()
% this function just for scatter all timeserises 
DU = 2.6870e+16;
save_fig = 0;
%load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2.mat');
%load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews_2017_high_quality.mat');% this dataset extended to 2017
%load('C:\Users\ZhaoX\Documents\paper\CF\fig\gbs_saoz_brewer_merra2_ews_2017_high_quality_dmp.mat');
load('C:\Users\ZhaoX\Documents\paper\CF\fig\gbs_saoz_brewer_merra2_ews_2017_high_quality_dmp_nogbs2017spring.mat');

SAOZ_V3.UTC = datenum(datetime(SAOZ_V3.Year,SAOZ_V3.Month,SAOZ_V3.Day));

%% extra QC for GBS
[GBS_CF,N_all,N_removed] = remove_single_measurements(GBS_CF);
disp(['GBS_C_F has ' num2str(N_all) ' measurements in total']);
disp(['GBS_C_F has ' num2str(N_removed) ' measurements are single and removed']);
GBS_CF = table2timetable(GBS_CF,'RowTimes','UTC_str');
GBS_CF.year = double(string(GBS_CF.year));
GBS_CF_hourly = retime(GBS_CF,'hourly','mean');
TF = isnan(GBS_CF_hourly.mean_vcd);
GBS_CF_hourly(TF,:)=[];


[GBS,N_all,N_removed] = remove_single_measurements(GBS);
disp(['GBS has ' num2str(N_all) ' measurements in total']);
disp(['GBS has ' num2str(N_removed) ' measurements are single and removed']);
GBS = table2timetable(GBS,'RowTimes','UTC_str');
GBS.year = double(string(GBS.year));
GBS_hourly = retime(GBS,'hourly','mean');
TF = isnan(GBS_hourly.mean_vcd);
GBS_hourly(TF,:)=[];


C = innerjoin(GBS_CF_hourly,GBS_hourly);
x = C.mean_vcd_GBS_CF_hourly./DU;
y = C.mean_vcd_GBS_hourly./DU;
line_fits(x,y);
dscatter(x,y);
colormap;
xlabel(['GBS_C_F TCO [DU]']);
ylabel(['GBS TCO [DU]']);

%% extra QC for SAOZ
[SAOZ_CF,N_all,N_removed] = remove_single_measurements(SAOZ_CF);
disp(['SAOZ_C_F has ' num2str(N_all) ' measurements in total']);
disp(['SAOZ_C_F has ' num2str(N_removed) ' measurements are single and removed']);
SAOZ_CF = table2timetable(SAOZ_CF,'RowTimes','UTC_str');
SAOZ_CF.year = double(string(SAOZ_CF.year));
SAOZ_CF_hourly = retime(SAOZ_CF,'hourly','mean');
TF = isnan(SAOZ_CF_hourly.mean_vcd);
SAOZ_CF_hourly(TF,:)=[];


[SAOZ,N_all,N_removed] = remove_single_measurements(SAOZ);
disp(['SAOZ has ' num2str(N_all) ' measurements in total']);
disp(['SAOZ has ' num2str(N_removed) ' measurements are single and removed']);
SAOZ = table2timetable(SAOZ,'RowTimes','UTC_str');
SAOZ.year = double(string(SAOZ.year));
SAOZ_hourly = retime(SAOZ,'hourly','mean');
TF = isnan(SAOZ_hourly.mean_vcd);
SAOZ_hourly(TF,:)=[];


C = innerjoin(SAOZ_CF_hourly,SAOZ_hourly);
x = C.mean_vcd_SAOZ_CF_hourly./DU;
y = C.mean_vcd_SAOZ_hourly./DU;
line_fits(x,y);
dscatter(x,y);
colormap;
xlabel(['SAOZ_C_F TCO [DU]']);
ylabel(['SAOZ TCO [DU]']);
%% fig 1
figure; hold all;
% MERRA-2
m = scatter(MERRA2.UTC,MERRA2.MERRA_Ozone,6,[0.5,0.5,0.5],'s','filled');
alpha(m,0.5);
% BrewerDS
TF = strcmp(Brewer.ObsCode,'DS');
b =scatter(Brewer.UTC(TF,:),Brewer.ColumnO3(TF,:),6,'s','filled');
alpha(b,0.5);
% GBS and GBS-CF
scatter(GBS.UTC,GBS.mean_vcd./DU,4,'o','filled');
alpha(0.5);
scatter(GBS_CF.UTC,GBS_CF.mean_vcd./DU,4,'*');
alpha(0.5);
% SAOZ, SAOZ-CF, and SAOZ-V3
scatter(SAOZ.UTC,SAOZ.mean_vcd./DU,4,'o','filled');
alpha(0.5);
scatter(SAOZ_CF.UTC,SAOZ_CF.mean_vcd./DU,4,'*');
alpha(0.5);
scatter(SAOZ_V3.UTC,SAOZ_V3.O3sr,4,'x');
alpha(0.5);
scatter(SAOZ_V3.UTC,SAOZ_V3.O3ss,4,'x');
alpha(0.5);

xlim([734130 737070]);
ylim([200 600]);
datetick('x','yyyy-mm','keeplimits');
ylabel('TCO [DU]');
legend('MERRA-2', 'Brewer_D_S', 'GBS','GBS_C_F','SAOZ','SAOZ_C_F','SAOZ_V_3');
grid on;
print_setting(1/2,save_fig,['Time_series_all']);

% figure 2
figure; hold all;

h1 = subplot(3,1,1);hold all;
% MERRA-2
m = scatter(MERRA2.UTC,MERRA2.MERRA_Ozone,6,[0.5,0.5,0.5],'s','filled');
alpha(m,0.5);
% BrewerDS
TF = strcmp(Brewer.ObsCode,'DS');
b =scatter(Brewer.UTC(TF,:),Brewer.ColumnO3(TF,:),6,'s','filled');
alpha(b,0.5);
xlim([734130 737070]);
ylim([200 600]);
ylabel('TCO [DU]');
grid on;

h2 = subplot(3,1,2);hold all;
% MERRA-2
m = scatter(MERRA2.UTC,MERRA2.MERRA_Ozone,6,[0.5,0.5,0.5],'s','filled');
alpha(m,0.5);
% GBS and GBS-CF
scatter(GBS.UTC,GBS.mean_vcd./DU,4,'o','filled');
alpha(0.5);
scatter(GBS_CF.UTC,GBS_CF.mean_vcd./DU,4,'*');
alpha(0.5);
xlim([734130 737070]);
ylim([200 600]);
ylabel('TCO [DU]');
grid on;


h3 = subplot(3,1,3);hold all;
% MERRA-2
m = scatter(MERRA2.UTC,MERRA2.MERRA_Ozone,6,[0.5,0.5,0.5],'s','filled');
alpha(m,0.5);
% SAOZ, SAOZ-CF
scatter(SAOZ.UTC,SAOZ.mean_vcd./DU,4,'o','filled');
alpha(0.5);
scatter(SAOZ_CF.UTC,SAOZ_CF.mean_vcd./DU,4,'*');
alpha(0.5);
xlim([734130 737070]);
ylim([200 600]);
ylabel('TCO [DU]');
grid on;

datetick('x','yyyy-mm','keeplimits');
print_setting(1/2,save_fig,['Time_series_all_subplots']);

function [data,N_all,N_removed] = remove_single_measurements(data)

data.daynum = fix(data.UTC);
TF_singleday = false(numel(data.daynum),1);
for i =1:numel(data.daynum)
    TF = data.daynum(i) == data.daynum;
    if sum(TF) == 1;
        TF_singleday(i) = true;
    end
end
N_all = height(data);
N_removed = sum(TF_singleday);
data(TF_singleday,:)=[];
