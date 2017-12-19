function general_output = X_days_low_f_test()
DU = 2.6870e+16;
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2_cf\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2_brewerzs\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_v2_cf_brewerzs\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_clear_test\'];
%datapath = ['E:\H\work\Eureka\GBS\CI\weather_impact_clear_cf_compare\'];
datapath = ['E:\H\work\Eureka\GBS\CI\archive\uncertainty\'];

%datapath = ['E:\H\work\Eureka\SAOZ\weather_impact_v2\'];
%datapath = ['E:\H\work\Eureka\SAOZ\weather_impact_v2_cf\'];
%datapath = ['E:\H\work\Eureka\SAOZ\weather_impact_v2_clear_test2\'];
%load([datapath 'weather_impact.mat']);
%load([datapath 'compare.mat']);
%load('E:\H\work\Eureka\GBS\CI\weather_impact_v2\GBS_Brewer_MERRA2');
%load('GBS_VCD_2010_2017_MERRA2_2010_2015_EWS.mat');
load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews_2017_high_quality.mat');


weather_types = {'All'};
%weather_types = {'Clear & Mainly Clear'};
% make groups of instrument/model pairs we are interested
instrument_pairs = {'GBS_Brewer_EWS','GBS_CF_Brewer_EWS','SAOZ_Brewer_EWS','SAOZ_CF_Brewer_EWS', ...
    'SAOZ_V3_reformat_Brewer_EWS','GBS_MERRA2_EWS','GBS_CF_MERRA2_EWS','SAOZ_MERRA2_EWS',...
    'SAOZ_CF_MERRA2_EWS','SAOZ_V3_reformat_MERRA2_EWS'};
%instrument_pairs = {'SAOZ_V3_reformat_MERRA2_EWS'};
%instrument_pairs = {'SAOZ_V3_reformat_Brewer_EWS'};
uncertainties = table;% this is a detailed uncertainty table, which has calculated uncertainties with using differen X-days of mean as low-f
general_output = table;% this is a simplified uncertainty table, which only has calculated uncertainties when X-days = 5 days
for i = 1:numel(instrument_pairs)
    for j = 1:numel(weather_types)
        
        plot_path = [datapath cell2mat(instrument_pairs(i)) '\' cell2mat(weather_types(j)) '\']; % prepare plot path
        mkdir(plot_path);
        cd(plot_path);
        
        uncertainties = [];
        eval(['gbs_X =' cell2mat(instrument_pairs(i)) ';']); % assign a table to gbs_X, which will be used in uncertianty estimation function
        
        if ~strcmp('All',weather_types(j)) % filter data by weather type
            weathers = split(weather_types(j),' & ');
            TF = strcmp(gbs_X.weather_median_ampm,cell2mat(weathers(1))) | strcmp(gbs_X.weather_median_ampm,cell2mat(weathers(2)));
            gbs_X(~TF,:) = [];
        end
        
        if strfind(cell2mat(instrument_pairs(i)),'MERRA2')>0 %
            gbs_X.mean_ColumnO3 = gbs_X.MERRA2_Ozone;% replace Brewer ozone with MERRA-2 ozone
        end
        %gbs_X.mean_vcd = gbs_X.MERRA2_Ozone.*DU;% replace GBS ozone with MERRA-2 ozone
        %TF = (strcmp(gbs_X.weather_median_ampm,'Cloudy')) ;
        %gbs_X(~TF,:) = [];
        
        gbs_vcd_type = 'normal';
        save_fig = 1;
        fig_title = [replace(cell2mat(instrument_pairs(i)),'_','-')];
        for step = 7
            try
                uncertainties_xstep = uncertainty_v4(gbs_X,gbs_vcd_type,0,fig_title,step);
                if step == 1
                    uncertainties = uncertainties_xstep;
                else
                    uncertainties = [uncertainties;uncertainties_xstep];
                end
            catch
                disp(['uncertainty estimation failled for ' cell2mat(instrument_pairs(i)) ' ' cell2mat(weather_types(j))]);
            end
            
        end
        
        steps = 7;
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
        legends = split(cell2mat(instrument_pairs(i)),'_');
        legend(cell2mat(legends(1)),cell2mat(legends(end-1)),'R value');
        title([cell2mat(weather_types(j))]);
        %print_setting('narrow2',save_fig,['uncertainty_estimation_Xdays_' cell2mat(instrument_pairs(i))]);
        print_setting('narrow2',0,['uncertainty_estimation_Xdays_' cell2mat(instrument_pairs(i))]);
        
        save('uncertainty_estimation');
        close all;
        
        %if ~strcmp('All',weather_types(j)) % filter data by weather type
            general_output(i,:) = uncertainties(1,:);
            general_output.Properties.RowNames(i) = instrument_pairs(i);
            
        %end
    end
end
