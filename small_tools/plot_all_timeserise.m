function scatter_all_timeserise()
% this function just for scatter all timeserises 
DU = 2.6870e+16;
save_fig = 1;
load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2.mat');

SAOZ_V3.UTC = datenum(datetime(SAOZ_V3.Year,SAOZ_V3.Month,SAOZ_V3.Day));

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
