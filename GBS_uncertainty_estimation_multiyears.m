function GBS_uncertainty_estimation_multiyears()

%year_list = [2010,2011,2013,2014,2015,2016];
year_list = [2011,2013,2014,2015,2016];
N_year = size(year_list);
working_dir = 'H:\work\Eureka\GBS\CI\2010_2016_uncertainty\';
cd(working_dir);
save_fig = 1;
%% get multiyear dataset
for i =1:1:N_year(2)
    year = num2str(year_list(i));
    data_paths = ['H:\work\Eureka\GBS\CI\' year '\CF_450_550_minCI\vs_Brewer\'];
    
    load([data_paths 'gbs_brewer.mat']);
        
    
    if i == 1
        gbs_brewer_multiyears = gbs_brewer;
        gbs_brewerzs_multiyears = gbs_brewerzs;
        
        gbscf_brewer_multiyears = gbscf_brewer;
        gbscf_brewerzs_multiyears = gbscf_brewerzs;
        
        list_HQ_day_multiyears = list_HQ_day;
    else
        gbs_brewer_multiyears = [gbs_brewer_multiyears;gbs_brewer];
        gbs_brewerzs_multiyears = [gbs_brewerzs_multiyears;gbs_brewerzs];
        
        gbscf_brewer_multiyears = [gbscf_brewer_multiyears;gbscf_brewer];
        gbscf_brewerzs_multiyears = [gbscf_brewerzs_multiyears;gbscf_brewerzs];
        
        list_HQ_day_multiyears = [list_HQ_day_multiyears;list_HQ_day];
    end   
end

clearvars -except  save_fig gbs_brewer_multiyears gbs_brewerzs_multiyears gbscf_brewer_multiyears gbscf_brewerzs_multiyears list_HQ_day_multiyears;

%% perform QC for GBS data then save
std_vcd_threshold_DS = mean(gbs_brewer_multiyears.std_vcd) + 3*std(gbs_brewer_multiyears.std_vcd);
mean_StdDevO3_threshold_DS = mean(gbs_brewer_multiyears.mean_StdDevO3) + 3*std(gbs_brewer_multiyears.mean_StdDevO3);
std_vcd_threshold_ZS = mean(gbs_brewerzs_multiyears.std_vcd) + 3*std(gbs_brewerzs_multiyears.std_vcd);
mean_StdDevO3_threshold_ZS = mean(gbs_brewerzs_multiyears.mean_StdDevO3) + 3*std(gbs_brewerzs_multiyears.mean_StdDevO3);

gbs_brewer_multiyears = QC_GBS_VCD(gbs_brewer_multiyears,std_vcd_threshold_DS,mean_StdDevO3_threshold_DS);
gbs_brewerzs_multiyears = QC_GBS_VCD(gbs_brewerzs_multiyears,std_vcd_threshold_ZS,mean_StdDevO3_threshold_ZS);
gbscf_brewer_multiyears = QC_GBS_VCD(gbscf_brewer_multiyears,std_vcd_threshold_DS,mean_StdDevO3_threshold_DS);
gbscf_brewerzs_multiyears = QC_GBS_VCD(gbscf_brewerzs_multiyears,std_vcd_threshold_ZS,mean_StdDevO3_threshold_ZS);

save('GBS_vs_Brewer.mat');


%% perform uncertainty estimation

fig_title = 'GBS vs BrewerDS';
uncertainties_gbs_brewer = uncertainty_v2(gbs_brewer_multiyears,save_fig,fig_title);

fig_title = 'GBS vs BrewerZS';
uncertainties_gbs_brewerzs = uncertainty_v2(gbs_brewerzs_multiyears,save_fig,fig_title);

fig_title = 'GBS-CF vs BrewerDS';
uncertainties_gbscf_brewer = uncertainty_v2(gbscf_brewer_multiyears,save_fig,fig_title);

fig_title = 'GBS-CF vs BrewerZS';
uncertainties_gbscf_brewerzs = uncertainty_v2(gbscf_brewerzs_multiyears,save_fig,fig_title);

save('GBS_vs_Brewer.mat');

%%
function data = QC_GBS_VCD(data,std_vcd_threshold,mean_StdDevO3_threshold)
N_initial = size(data);
%% for GBS: we will filter out any data has std_vcd larger than mean + 3 sigma value, typically will filter out 1% of data
TF1 = data.std_vcd > std_vcd_threshold;

%% for Brewer: we will filter out any data has std_vcd larger than mean + 3 sigma value, typically will filter out 1% of data
TF2 = data.mean_StdDevO3 > mean_StdDevO3_threshold;

%% for GBS: we only want measurements that have both sr and ss measurements
%TF2 = isnan(data.sigma_mean_vcd); 
j = 0; TF4(1) = false;
day_stamp = fix(data.LTC);
for i =1:1:N_initial(1)
    TF3 = day_stamp(i) == day_stamp;
    if sum(TF3) == 1
        %single_date = datestr(data.LTC);
        %disp(['find a single measuremnt at ' ]);
        %disp(single_date);
        j = j + 1;
        TF4(i) = true;
    elseif sum(TF3) == 2
        TF4(i) = false;
    end
end

TF = TF1 | TF2 | TF4';

data(TF,:) = [];

p_filtered = sum(TF)/N_initial(1)*100;
disp([num2str(p_filtered) '% data been filtered in QC ...']);




