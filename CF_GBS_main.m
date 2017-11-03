function CF_GBS_main()
% this is the main function to
% 1) read in GBS QDOAS output, then
% 2) lable the data with sky flags
% 3) the labled data will be processed by GBS VCD
% package, which will give both cloud filtered and normal ozone VCDs.
% 4) the GBS VCDs will be paried with Brewer data (both ZS and DS) and
% perform random uncertainty anaylsis

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

% read the input txt file, please make sure you have a pareller folder to
% this 'CF_package', which is named as 'CF_package_local'! This local
% folder will content input files and other files that we won't share on
% Github!
input_table = read_input_file();

% adding matlab searching path
addpath(input_table.CF_package_path)
addpath(input_table.Supporting_code_path);
addpath(input_table.CF_package_local_file_path);

%% interpret input.table to simple input values %%
% year label for data will be processed
year = input_table.year;

% read L2 data
% CF package only can read preprocessed L2 mat file!
data_input =  input_table.data_input;
load(data_input);

% need check the unit for O4 in the dSCDs file(eg, GBS v2 data)
% o4_xs_unit_convertion = 1;% 1= *1e40; % if use Greenblatt o4.xs, need *1e40 to convert unit to [cm5/molec]
o4_xs_unit_convertion = input_table.o4_xs_unit_convertion;

% load sonde data
load(input_table.data_sonde);

% EWS weather record
EWS_raw_data_folder = strcat(input_table.EWS_folder,year); 

% load Brewer data
load(input_table.Brewer_data_path); 
eval(['brewer_all = combined_raw_' year ';']);


% add VCD code package
VCD_code_path = input_table.VCD_code_path;
addpath(VCD_code_path);

% this is add-on, which will calculate time of the ref spec.
addpath([ VCD_code_path 'time_from_SZA']);

% define plot path
plot_path = input_table.plot_path;

% auto save figure or not: 0 = no; 1 = yes
save_fig = input_table.save_fig;

% if you want escape some steps, in the input file, set escape to true! use with caution!
if input_table.escape == true
    load([plot_path 'temp.mat']);
    input_table.escape = true;
    try
        disp('temp file loaded ... ');
        disp(['CF package will start from Step: ' num2str(step_number + 1)]);
    catch
        disp('Warning: temp file loading error/step number missing!');
        disp('CF package will make a fresh run this time');
        step_number = 0;
        input_table.escape = false;
    end
else
    disp('CF package will make a fresh run');
    step_number = 0;
end

%% 1) read in QDOAS data
if (input_table.escape == false) | (step_number < 1) % check if we want escape this step
    try       
        if strcmp(input_table.instrument,'GBS')
            % for UT-GBS
            data = ci_o4_testing_plot_v2(data_input,year,plot_path,save_fig);
            if o4_xs_unit_convertion == 1
                data.O3SlColo4 = data.O3SlColo4.*1e40;
                data.O3SlErro4 = data.O3SlErro4.*1e40;
            end
        elseif strcmp(input_table.instrument,'SAOZ')
            %for SAOZ
            eval(['data = saoz' year ';']);
        end
        
        disp('>>> Step 1 finished');
        step_number = 1;
        cd(plot_path);
        save('temp.mat');
        
    catch
        disp('Warning: step 1 failed');
        cd(plot_path);
        save('temp.mat');
        email_notice('CODE STOPPED @ STEP-1 DUE TO ERROR');
        pause;
    end
end

%% 2) sky condtion
if (input_table.escape == false) | (step_number < 2)
    try
        %% 2.1 ) read in EWS data  
        EWS = make_EWS_weather_table_v2(EWS_raw_data_folder);
        data = pair_CI_weather_v2(data,EWS);

        %% 2.2) use the cloudy bin method to label sky condition
        % add calibrated CI
        data = fit_cloudy_bin_v3_450_550(data,plot_path,save_fig);
        % add TSI 
        [data_raw,data,list_HQ_day] = smooth_flag_v2(data,save_fig);
        save('dSCDs_CF.mat','data');
        
        disp('>>> Step 2 finished');
        step_number = 2;

    catch
        disp('Warning: step 2 failed');
        cd(plot_path);
        save('temp.mat');
        email_notice('CODE STOPPED @ STEP-2 DUE TO ERROR');
        pause;
    end
    save('temp.mat');
    disp('Step 2 finished ... wait for 60 s ... ');pause(60);
end

%% 3) process sky flaged dSCDs in our modified VCD package
if (input_table.escape == false) | (step_number < 3)
    try
        year_num = str2num(year);
        QDOAS_data_dir = plot_path;
        QDOAS_data_file = 'dSCDs_CF.mat';
        
        % use VCD package to convert dSCDs to VCD
        [VCD,VCD_CF, dscd_S, dscd_SCF, rcd_S, rcd_SCF, avg_vcd, avg_vcdCF, qdoas_filt] = DSCD_to_VCD(year,VCD_code_path,plot_path,save_fig,QDOAS_data_dir,QDOAS_data_file,sonde);
        
        cd(plot_path);
        % plot timeserise (brewer DS, ZS, and the instrument using this CF package)
        plot_timeserise(VCD,VCD_CF,brewer_all,plot_path,save_fig);
        save('temp.mat');
        
        disp('>>> Step 3 finished');
        step_number = 3;
        
    catch
        disp('Warning: step 3 failed');
        cd(plot_path);
        save('temp.mat');
        email_notice('CODE STOPPED @ STEP-3 DUE TO ERROR');
        pause;
    end
 
end
%% 4) vs Brewer
if (input_table.escape == false) | (step_number < 4)
    try
        % pair Brewer DS data with GBS/SAOZ data
        fig_name = '';
        ds_zs = 0;
        gbs_brewer = pair_Brewer_GBS_daily_v2(brewer_all,ds_zs,VCD,save_fig,fig_name,plot_path);
        
        % pair Brewer DS data with GBS/SAOZ cloud filtered data
        fig_name = '-CF';
        ds_zs = 0;
        gbscf_brewer = pair_Brewer_GBS_daily_v2(brewer_all,ds_zs,VCD_CF,save_fig,fig_name,plot_path);

        % pair Brewer ZS data with GBS/SAOZ data
        fig_name = '';
        ds_zs = 1;
        gbs_brewerzs = pair_Brewer_GBS_daily_v2(brewer_all,ds_zs,VCD,save_fig,fig_name,plot_path);

        % pair Brewer ZS data with GBS/SAOZ cloud filtered data
        fig_name = '-CF';
        ds_zs = 1;
        gbscf_brewerzs = pair_Brewer_GBS_daily_v2(brewer_all,ds_zs,VCD_CF,save_fig,fig_name,plot_path);
        
        disp('>>> Step 4 finished');
        step_number = 4;
        
        cd(plot_path);
        save('temp.mat');

    catch
        disp('Warning: step 4 failed');
        cd(plot_path);
        save('temp.mat');
        email_notice('CODE STOPPED @ STEP-4 DUE TO ERROR');
        pause;
    end
end


%% 5) uncertaity estimation
if (input_table.escape == false) | (step_number < 5)
    try
        % calculate uncertainties
        fig_title = 'BrewerDS vs GBS';
        gbs_vcd_type = 'normal';
        uncertainties_gbs_brewer = uncertainty_v3(gbs_brewer,gbs_vcd_type,save_fig,fig_title);

        fig_title = 'BrewerDS vs GBS(langley)';
        gbs_vcd_type = 'langley';
        uncertainties_gbsl_brewer = uncertainty_v3(gbs_brewer,gbs_vcd_type,save_fig,fig_title);

        fig_title = 'BrewerDS vs GBS-CF';
        gbs_vcd_type = 'normal';
        uncertainties_gbscf_brewer = uncertainty_v3(gbscf_brewer,gbs_vcd_type,save_fig,fig_title);

        fig_title = 'BrewerDS vs GBS-CF(langley)';
        gbs_vcd_type = 'langley';
        uncertainties_gbscfl_brewer = uncertainty_v3(gbscf_brewer,gbs_vcd_type,save_fig,fig_title);

        fig_title = 'BrewerZS vs GBS';
        gbs_vcd_type = 'normal';
        uncertainties_gbs_brewerzs = uncertainty_v3(gbs_brewerzs,gbs_vcd_type,save_fig,fig_title);

        fig_title = 'BrewerZS vs GBS(langley)';
        gbs_vcd_type = 'langley';
        uncertainties_gbsl_brewerzs = uncertainty_v3(gbs_brewerzs,gbs_vcd_type,save_fig,fig_title);

         fig_title = 'BrewerZS vs GBS-CF';
        gbs_vcd_type = 'normal';
        uncertainties_gbscf_brewerzs = uncertainty_v3(gbscf_brewerzs,gbs_vcd_type,save_fig,fig_title);

        fig_title = 'BrewerZS vs GBS-CF(langley)';
        gbs_vcd_type = 'langley';
        uncertainties_gbscfl_brewerzs = uncertainty_v3(gbscf_brewerzs,gbs_vcd_type,save_fig,fig_title);

        uncertainties = [uncertainties_gbs_brewer;uncertainties_gbscf_brewer;uncertainties_gbs_brewerzs;uncertainties_gbscf_brewerzs;uncertainties_gbsl_brewer;uncertainties_gbscfl_brewer;uncertainties_gbsl_brewerzs;uncertainties_gbscfl_brewerzs];
        
        disp('>>> Step 5 finished');
        step_number = 5;
        cd(plot_path);
        save('temp.mat');
        
    catch
        disp('Warning: step 5 failed');
        cd(plot_path);
        save('temp.mat');
        email_notice('CODE STOPPED @ STEP-5 DUE TO ERROR');
        pause;
    end

end

%% 6) save data
gbs = data;
brewer = brewer_all;

save('gbs_brewer.mat','gbs','brewer','EWS','gbs_brewer','gbscf_brewer','gbs_brewerzs','gbscf_brewerzs','VCD','VCD_CF','list_HQ_day','qdoas_filt','rcd_S','rcd_SCF','year','uncertainties');
try
    email_notice('CODE SUCCESSFULLY FINISHED!');
catch
    disp('CAN NOT SEND EMAIL, PLS CHECK SETTINGS IN "email_notice".');
end
close all;
