function gbs_brewer = pair_Brewer_GBS_daily_v2(brewer_all,ds_zs,VCD,save_fig,fig_name,plot_path)
size = 1;
%save_fig = 1;
%fig_name = ['BrewerZS vs GBS'];
%plot_path = 'H:\work\Eureka\GBS\CI\2011\UTGBS\vs_Brewer';
plot_path = [plot_path 'vs_Brewer'];
mkdir(plot_path);
cd(plot_path);
DU = 2.69e16;
%% load Brewer
%load('H:\work\Eureka\Brewer\Bruker69_2004_2014_all_modes.mat');
%brewer_all = combined_raw_2011;
brewer_all.LTC = brewer_all.UTC - 5/24;
TF = strcmp(brewer_all.ObsCode, 'DS');
brewer_ds = brewer_all(TF,:);
TF = strcmp(brewer_all.ObsCode, 'ZS');
brewer_zs = brewer_all(TF,:);
if ds_zs == 0
    brewer = brewer_ds; % use brewer direct-sun data
    fig_name = ['BrewerDS vs GBS' fig_name];
elseif ds_zs == 1
    brewer = brewer_zs; % use brewer zenith-sky data
    fig_name = ['BrewerZS vs GBS' fig_name];
end
try
    brewer.Time = [];
catch
end
brewer.ObsCode = [];

%% load GBS
%load('H:\work\Eureka\GBS\CI\2011\UTGBS\VCD3\matlab_2.mat');
gbs = VCD;
%gbs = VCD_table; % unfiltered GBS data
%gbs = VCD_table2; % filtered GBS data
gbs.UTC = fd_to_UTC(gbs.year,gbs.fd);
gbs.LTC = gbs.UTC - 5/24;


%% average Brewer to am/pm and pair with GBS
%N = size(gbs);
N = height(gbs);
gbs_matched = logical(zeros(N(1),1));
j = 1; 
%ampm_index = 0;
for i = 1:1:N(1)
    %TF = abs(brewer.LTC - gbs.LTC(i)) < 0.5;
    if gbs.ampm(i) == 0 % for moring measurements
        start_time = fix(gbs.LTC(i));
    else % for afternoon measurements
        start_time = fix(gbs.LTC(i))+0.5;    
    end
    %disp(['start time: ' datestr(start_time)]);
    end_time = start_time + 0.5;
    %disp(['end time: ' datestr(end_time)]);
    TF = (brewer.LTC > start_time) & (brewer.LTC < end_time);
    if sum(TF) > 0
        func = @mean;
        brewer_avg(j,:) = varfun(func,brewer(TF,:));
        %ampm_index(j,:) = gbs.ampm(i);
        gbs_matched(i,:) = true;
        j = j + 1;      
    end
end

if sum(gbs_matched) == 0
    disp('Warning: No Brewer measurements found coincident with GBS measurements!');
    disp('Escape from this step ... ');
    gbs_brewer = 0;
else
    gbs_brewer = [gbs(gbs_matched,:) , brewer_avg];
    
    
    TF = isnan(gbs_brewer.mean_vcd);
    gbs_brewer(TF,:) = [];
    N_gbs_brewer = height(gbs_brewer);
    %% plots
    x = gbs_brewer.mean_vcd./DU; % GBS ozone
    x_err = mean(gbs_brewer.std_vcd)./DU;
    p_x_err = x_err/mean(x)*100;
    y = gbs_brewer.mean_ColumnO3; % Brewer ozone
    y_err = mean(gbs_brewer.mean_StdDevO3);
    p_y_err = y_err/mean(y)*100;
    
    linear_fits(x,y,x_err,y_err);
    print_setting(size,0,fig_name);
    textbp(['N = ' num2str(N_gbs_brewer(1))]);
    %textbp(['GBS_e_r_r(std) = ' num2str(mean(x_err)) '(' num2str(p_x_err)  '%)' ' ; ' 'Brewer_e_r_r(std) = ' num2str(mean(y_err))  '(' num2str(p_y_err)  '%)' ]);
    textbp(['GBS_e_r_r(std) = ' num2str(mean(x_err)) '(' num2str(p_x_err)  '%)' ]);
    textbp(['Brewer_e_r_r(std) = ' num2str(mean(y_err))  '(' num2str(p_y_err)  '%)' ]);
    title(fig_name);
    xlabel('GBS TCO [DU]');
    ylabel('Brewer TCO [DU]');
    print_setting(size,save_fig,fig_name);
    
end

