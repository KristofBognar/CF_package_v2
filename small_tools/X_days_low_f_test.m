function uncertainties = X_days_low_f_test()
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2_cf\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2_brewerzs\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2_cf_brewerzs\'];

%datapath = ['E:\H\work\Eureka\SAOZ\weather_impact_v2\'];
datapath = ['E:\H\work\Eureka\SAOZ\weather_impact_v2_cf\'];
load([datapath 'weather_impact.mat']);
cd(datapath);

uncertainties = table;
gbs_brewer = final_table_concat;
gbs_vcd_type = 'normal';
save_fig = 1;
%fig_title = ['GBS(2010-2017)'];
fig_title = ['SAOZ(2010-2017)'];
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
legend('SAOZ-CF','Brewer DS','R value');
print_setting('narrow2',save_fig,['uncertainty_estimation_Xdays_temp']);


save('uncertainty_estimation');
