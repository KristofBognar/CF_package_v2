function uncertainties = X_days_low_f_test()
DU = 2.6870e+16;
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2_cf\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2_brewerzs\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2_cf_brewerzs\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_clear_test\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_clear_cf_compare\'];
datapath = ['E:\H\work\Eureka\GBS\CI\MERRA2\'];

%datapath = ['E:\H\work\Eureka\SAOZ\weather_impact_v2\'];
%datapath = ['E:\H\work\Eureka\SAOZ\weather_impact_v2_cf\'];
%datapath = ['E:\H\work\Eureka\SAOZ\weather_impact_v2_clear_test2\'];
%load([datapath 'weather_impact.mat']);
%load([datapath 'compare.mat']);
%load('E:\H\work\Eureka\GBS\CI\weather_impact_v2\GBS_Brewer_MERRA2');
load('GBS_VCD_2010_2017_MERRA2_2010_2015_EWS.mat');
cd(datapath);

uncertainties = table;
%gbs_brewer = final_table_concat;
%gbs_brewer = data_cfonly;
gbs_brewer = GBS_VCD_MERRA2_EWS;
gbs_brewer(isnan(gbs_brewer.MERRA2_Ozone),:) = [];
gbs_brewer.mean_ColumnO3 = gbs_brewer.MERRA2_Ozone;% replace Brewer ozone with MERRA-2 ozone
%gbs_brewer.mean_vcd = gbs_brewer.MERRA2_Ozone.*DU;% replace GBS ozone with MERRA-2 ozone
TF = (strcmp(gbs_brewer.weather_median_ampm,'Cloudy')) ;
gbs_brewer(~TF,:) = [];

gbs_vcd_type = 'normal';
save_fig = 1;
fig_title = ['GBS(2010-2015) GBS(Cloudy) vs MERRA2'];
%fig_title = ['SAOZ(2010-2017)'];
for step = 1:1:10
    uncertainties_xstep = uncertainty_v4(gbs_brewer,gbs_vcd_type,save_fig,fig_title,step);
    if step == 1
        uncertainties = uncertainties_xstep;
    else
        uncertainties = [uncertainties;uncertainties_xstep];
    end
end

steps = 1:1:10;
figure;hold all;
yyaxis left
errorbar(steps,uncertainties.pu_GBS,uncertainties.pu_GBS_err);
errorbar(steps,uncertainties.pu_Brewer,uncertainties.pu_Brewer_err);
ylabel(['TCO uncertainty [%]']);
ylim([0 5]);
yyaxis right
plot(steps,uncertainties.rho,'s');
ylabel('R');
ylim([0 1]);
xlabel('X days (low-f ozone average)');

%legend('GBS','Brewer DS','R value');
%legend('GBS','Brewer ZS','R value');
%legend('GBS-CF','Brewer ZS','R value');
%legend('SAOZ','Brewer DS','R value');
%legend('SAOZ-CF','Brewer DS','R value');
%legend('GBS','MERRA-2','R value');
legend('GBS','MERRA2','R value');
print_setting('narrow2',save_fig,['uncertainty_estimation_Xdays_GBS_VCD_vs_MERRA2_2010_2015']);


save('uncertainty_estimation');
