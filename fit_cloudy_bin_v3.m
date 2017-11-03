function data_output = fit_cloudy_bin_v3(data_input,plot_path,save_fig)

%% input file %%
SZA_filter = 1; % note, the cloudy bin method needs use SZA filter, so set this to 1 = yes, unless you want test some feature
% the output has two format, "data" is SZA filtered results, "data_output"
% is unfiltered results (for example, if SZA_filter = 1, then data_output give calibrated CI and classification flags with use SZA filter, but the output data hold resulst from all SZA. For those measurements "should be" filtered by SZA, their sky flags are set to 0).

use_min_COD = 1; % 1 = use the lowest CI simulated in RMT as high COD value, 0 = use a single simulation result as high COD CI

%save_fig = 0;
fig_size = 1/2;
%plot_path = 'H:\work\Eureka\GBS\CI\2011\UTGBS\plots\CI_callibration\'

%data_input = 'H:\work\Eureka\Eureka_weather_station\UTGBS_CI_2011_weather_HQ.mat';
%data_input = 'H:\work\Eureka\Eureka_weather_station\UTGBS_CI_2010_weather_HQ.mat';
%load(data_input);
data = data_input;

RTM_data_path = 'H:\work\Eureka\GBS\CI\matlab\CF_package\database\';

%% low COD condition
%data_path = 'H:\work\Eureka\GBS\CI\2010\UTGBS\plots\RTM\cloud_tau2_1to3km\';
%data_path = 'H:\work\Eureka\GBS\CI\2010\UTGBS\plots\RTM\cloud_tau1dot5_2to3km_albedo_0dot06\';
data_path = [RTM_data_path 'cloud_tau1dot5_2to3km_albedo_0dot06\'];
simulation_type = 'lowCOD';
model_lowCOD = SCIATRAN2_CI(data_path,simulation_type,plot_path);
%% high COD condition
%data_path = 'H:\work\Eureka\GBS\CI\2010\UTGBS\plots\RTM\cloud_tau12_1to3km\';
%data_path = 'H:\work\Eureka\GBS\CI\2010\UTGBS\plots\RTM\cloud_tau9_2to3km_albedo_0dot06\';
if use_min_COD == 1
    data_path = RTM_data_path;
    simulation_type = 'highCOD';
    [model, model_min] = SCIATRAN2_min_CI(data_path,simulation_type,plot_path);
    model_highCOD = model_min;
else
    data_path = [RTM_data_path 'cloud_tau9_2to3km_albedo_0dot06\'];
    simulation_type = 'highCOD';
    model_highCOD = SCIATRAN2_CI(data_path,simulation_type,plot_path);
end
%% clear condition
%data_path = 'H:\work\Eureka\GBS\CI\2010\UTGBS\plots\RTM\clear_sky\';
data_path = [RTM_data_path 'clear_sky\'];
simulation_type = 'clear_sly';
model_clear = SCIATRAN2_CI(data_path,simulation_type,plot_path);
%% lowAOD condition
%data_path = 'H:\work\Eureka\GBS\CI\2010\UTGBS\plots\RTM\aerosold_lowtran_23km_albedo_0dot06\';
data_path = [RTM_data_path 'aerosold_lowtran_23km_albedo_0dot06\'];
simulation_type = 'lowAOD';
model_lowAOD = SCIATRAN2_CI(data_path,simulation_type,plot_path);
%% verylowhAOD condition
%data_path = 'H:\work\Eureka\GBS\CI\2010\UTGBS\plots\RTM\aerosold_lowtran_50km_albedo_0dot06\';
data_path = [RTM_data_path 'aerosold_lowtran_50km_albedo_0dot06\'];
simulation_type = 'verylowAOD';
model_verylowAOD = SCIATRAN2_CI(data_path,simulation_type,plot_path);


%data = UTGBS_CI;
%data.CI = UTGBS_CI.CI440550;
try
    data.CI = data.CI440550; % define the main CI will be used for couldy bin
end
data_output = data;


%% SZA_filter 
if SZA_filter == 1
    TF_SZA = data.SZA <= 85;
    data(~TF_SZA,:) = [];
    if sum(TF_SZA) == 0
        disp(['Warning: no measurements were made for SZA <= 85!']);
    end
end


N = size(data);

CI_lowCOD = spline(model_lowCOD.SZA,model_lowCOD.CI,data.SZA);
CI_highCOD = spline(model_highCOD.SZA,model_highCOD.CI,data.SZA);
CI_verylowAOD = spline(model_verylowAOD.SZA,model_verylowAOD.CI,data.SZA);

beta = 0.7:0.005:1.2;
N_beta = size(beta);
for i = 1:1:N_beta(2)
    cal_CI = beta(i).*data.CI;
    meas_diff_CI_lowCOD = cal_CI - CI_lowCOD;
    meas_diff_CI_highCOD = cal_CI - CI_highCOD;
    meas_diff_verylowAOD = cal_CI - CI_verylowAOD;

    TF = (meas_diff_CI_lowCOD <= 0) & (meas_diff_CI_highCOD >= 0);% circled in the "cloudy-box"
    p(i) = sum(TF)/N(1);
    TF2 = meas_diff_verylowAOD >= 0;%defined as good condition
    p2(i) = sum(TF2)/N(1);
    TF3 = (meas_diff_verylowAOD < 0) & (meas_diff_CI_lowCOD > 0);%defined as mediocre
    p3(i) = sum(TF3)/N(1);
    TF4 = meas_diff_CI_lowCOD <= 0;%defined as cloudy
    p4(i) = sum(TF4)/N(1);
end

[M,I] = max(p);
beta_max = beta(I);
cal_CI = beta_max.*data.CI;

%% output data
data.cal_CI = cal_CI;
disp(['Estimated \beta = ' num2str(beta_max)]);
%% output the CI-sky flag
meas_diff_CI_lowCOD = cal_CI - CI_lowCOD;
meas_diff_CI_highCOD = cal_CI - CI_highCOD;
meas_diff_verylowAOD = cal_CI - CI_verylowAOD;
TF = (meas_diff_CI_lowCOD <= 0) & (meas_diff_CI_highCOD >= 0);% circled in the "cloudy-box"
data.cloudybox = TF;
TF2 = meas_diff_verylowAOD >= 0;%defined as clear condition
data.clear = TF2;
TF3 = (meas_diff_verylowAOD < 0) & (meas_diff_CI_lowCOD > 0);%defined as mediocre
data.mediocre = TF3;
TF4 = meas_diff_CI_lowCOD <= 0;%defined as cloudy
data.cloudy = TF4;

%% output data (note, data_output retained all measurements (all SZA)!)
N_raw = size(data_output);
data_output.cal_CI = beta_max.*data_output.CI;
data_output.cloudybox = logical(zeros(N_raw(1),1));
data_output.clear = logical(zeros(N_raw(1),1));
data_output.mediocre = logical(zeros(N_raw(1),1));
data_output.cloudy = logical(zeros(N_raw(1),1));

data_output.cloudybox(TF_SZA,:) = data.cloudybox;
data_output.clear(TF_SZA,:) = data.clear;
data_output.mediocre(TF_SZA,:) = data.mediocre;
data_output.cloudy(TF_SZA,:) = data.cloudy;


%% plots
if sum(TF_SZA) > 0
figure;hold all;
plot(beta,p2*100,'.-');
plot(beta,p3*100,'.-');
plot(beta,p4*100,'.-');
plot(beta,p*100,'.-');
xlabel('\beta');
ylabel('% of measurements');
textbp(['Estimated \beta = ' num2str(beta_max)]);
legend('clear','mediocre','cloudy','cloudy-box');
print_setting(fig_size, save_fig, ['Beta_estimation']);



%% CI cloudy screen SZA density plot
figure;hold all;
dscatter(data.SZA,data.CI);
dscatter(data.SZA,data.cal_CI);
plot(model_lowCOD.SZA,model_lowCOD.CI,'--');
plot(model_highCOD.SZA,model_highCOD.CI,'--');
plot(model_highCOD.SZA,model_clear.CI,'');
plot(model_lowAOD.SZA,model_lowAOD.CI,'.-');
plot(model_verylowAOD.SZA,model_verylowAOD.CI,'.-');
legend('CI','calibrated CI','low COD','high COD','AOD/COD 0','low AOD','verylow AOD');
xlabel('SZA');
ylabel('Colour Index');
print_setting(fig_size,save_fig,['CI cloudy screen SZA density_calibration']);

%% CI cloudy screen SZA density plot
figure;hold all;
dscatter(data.SZA,data.CI);
dscatter(data.SZA,data.cal_CI);
gscatter(data.SZA,data.cal_CI,data.Weather_simple_clearL1);
gscatter(data.SZA,data.cal_CI,data.Weather_simple_clearL2);
gscatter(data.SZA,data.cal_CI,data.Weather_simple_clearL3);
gscatter(data.SZA,data.cal_CI,data.Weather_simple_cloudyL1);
gscatter(data.SZA,data.cal_CI,data.Weather_simple_cloudyL2);
plot(model_lowCOD.SZA,model_lowCOD.CI,'--');
plot(model_highCOD.SZA,model_highCOD.CI,'--');
plot(model_highCOD.SZA,model_clear.CI,'');
plot(model_lowAOD.SZA,model_lowAOD.CI,'.-');
plot(model_verylowAOD.SZA,model_verylowAOD.CI,'.-');
legend('CI','calibrated CI','Weather clearL1','Weather clearL1','Weather clearL2','Weather clearL2','Weather clearL3','Weather clearL3','Weather cloudyL1','Weather cloudyL1','Weather cloudyL2','Weather cloudyL2','low COD','high COD','AOD/COD 0','low AOD','verylow AOD');
xlabel('SZA');
ylabel('Colour Index');
print_setting(fig_size,save_fig,['CI cloudy screen SZA density_calibration_EWS']);

figure;hold all;bin = 0.02;
h_hist = histogram(data.CI);
h_hist.Normalization = 'probability';
h_hist.BinWidth = bin;
h_hist2 = histogram(data.cal_CI);
h_hist2.Normalization = 'probability';
h_hist2.BinWidth = bin;
h_hist3 = histogram(data.cal_CI(TF));
h_hist3.Normalization = 'probability';
h_hist3.BinWidth = bin;
h_hist4 = histogram(data.cal_CI(TF2));
h_hist4.Normalization = 'probability';
h_hist4.BinWidth = bin;
h_hist5 = histogram(data.cal_CI(TF3));
h_hist5.Normalization = 'probability';
h_hist5.BinWidth = bin;
h_hist6 = histogram(data.cal_CI(TF4));
h_hist6.Normalization = 'probability';
h_hist6.BinWidth = bin;
legend('CI','calibrated CI','cloudy-box','clear','mediocre','cloudy');
xlabel('CI');
ylabel('f (normalized)');
print_setting(fig_size,save_fig,['hist_normalized']);


figure;hold all;bin = 0.02;
h_hist = histogram(data.CI);
%h_hist.Normalization = 'probability';
h_hist.BinWidth = bin;
h_hist2 = histogram(cal_CI);
%h_hist2.Normalization = 'probability';
h_hist2.BinWidth = bin;
h_hist3 = histogram(data.cal_CI(TF));
%h_hist3.Normalization = 'probability';
h_hist3.BinWidth = bin;
h_hist4 = histogram(data.cal_CI(TF2));
%h_hist4.Normalization = 'probability';
h_hist4.BinWidth = bin;
h_hist5 = histogram(data.cal_CI(TF3));
%h_hist5.Normalization = 'probability';
h_hist5.BinWidth = bin;
h_hist6 = histogram(data.cal_CI(TF4));
%h_hist6.Normalization = 'probability';
h_hist6.BinWidth = bin;
legend('CI','calibrated CI','cloudy-box','clear','mediocre','cloudy');
xlabel('CI');
ylabel('f');
print_setting(fig_size,save_fig,['hist']);

end

