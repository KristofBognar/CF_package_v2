function uncertainties = uncertainty_v4(gbs_brewer,gbs_vcd_type,save_fig,fig_title,X_Day_step)

%% input
%save_fig = 1;
%fig_file = 'H:\work\Eureka\GBS\CI\2011\UTGBS\vs_Brewer\'
%mkdir(fig_file);
%load('H:\work\Eureka\GBS\Validation_ACE_OSIRIS\formated_data\gbs_and_saoz_paried.mat');
%load('H:\work\Eureka\GBS\CI\2011\UTGBS\vs_Brewer\matlab.mat');
DU=2.69e16;

%X_Day_step = 9; % X day mean, this is for calculate the low frequency (X day averged) ozone value
data = gbs_brewer;
try
    TF = str2num(data.year) >= 2010;
catch 
    data.year = num2str(datetime(datevec(data.UTC)).Year);
    data.day = data.DoY;
    TF = str2num(data.year) >= 2010;
end
data = data(TF,:);
TF = isnan(data.mean_ColumnO3) | isnan(data.mean_vcd);
data(TF,:) = [];

if ~istable(data)
    disp('Can not estimate uncertainties: No Brewer measurements found conincident with GBS measurements!');
    %uncertainties = 0;
    uncertainties = table;
    uncertainties.u_GBS = -9999;
    uncertainties.pu_GBS = -9999;
    uncertainties.u_Brewer = -9999;
    uncertainties.pu_Brewer = -9999;
    uncertainties.u_X = -9999;
    uncertainties.pu_X = -9999;
else
        %% M1 and M2 subtract daily mean, to get only high frequency part
    for year = min(str2num(data.year)):1:max(str2num(data.year))
        for day = min(data.day):1:max(data.day)
            TF = (str2num(data.year) == year) & (data.day == day);
            if sum(TF) == 1
                data(TF,:) = []; % we only keep days have both am. pm. values
            elseif sum(TF) == 2
                data.daily_gbs_mean(TF,:) = mean(data.mean_vcd(TF,:));
                try
                    data.daily_gbs_mean_langley(TF,:) = mean(data.langley_vcd(TF,:));
                catch
                    %disp('No langley slop VCD was found');
                end
                data.daily_brewer_mean(TF,:) = mean(data.mean_ColumnO3(TF,:));
            end
        end   
    end
    
    for year = min(str2num(data.year)):1:max(str2num(data.year))
        for X_day = min(data.day):X_Day_step:max(data.day)
            TF = (str2num(data.year) == year) & (data.day>=X_day) & (data.day<X_day+X_Day_step);
            if sum(TF) == 1
            %if (sum(TF) < fix(X_Day_step/2 +2)) | sum(TF) == 1
                data(TF,:) = []; % we only keep days have both am. pm. values
            elseif sum(TF) > 1
            %elseif sum(TF) >= fix(X_Day_step/2 +2)
                data.Xday_gbs_mean(TF,:) = mean(data.mean_vcd(TF,:));
                try
                    data.Xday_gbs_mean_langley(TF,:) = mean(data.langley_vcd(TF,:));
                catch
                    %disp('No langley slop VCD was found');
                end
                data.Xday_brewer_mean(TF,:) = mean(data.mean_ColumnO3(TF,:));
            end
        end
    end
        
    %%
    if strcmp(gbs_vcd_type,'normal')
        if max(data.mean_vcd) > 1e18
            M1 = (data.mean_vcd - data.Xday_gbs_mean)./DU; % GBS TCO using our normal routine (RCD -> VCDs -> mean VCD)
        else
            M1 = (data.mean_vcd - data.Xday_gbs_mean); % GBS TCO using our normal routine (RCD -> VCDs -> mean VCD)
        end
    elseif strcmp(gbs_vcd_type,'langley')
        M1 = (data.langley_vcd - data.Xday_gbs_mean_langley)./DU; % GBS TCO using slope of langley fit as VCD
    end
    %M2 = data.mean_ColumnO3 - data.daily_brewer_mean; % brewer TCO
    M2 = data.mean_ColumnO3 - data.Xday_brewer_mean; % brewer TCO

    %fig_title = 'Brewer vs GBS';
    
    %% filters: 
    %filter 1 : we will filter any measuremnts with abs(TCO - daily_mean) > 100 DU
    TF1 = abs(M1) > 100;
    TF2 = abs(M2) > 100;
    TF = TF1 | TF2;
    M1(TF,:) = [];
    M2(TF,:) = [];
    %filter 2 : we will filter TCO in table as NaN
    TF1 = isnan(M1);
    TF2 = isnan(M2);
    TF = TF1 | TF2;
    M1(TF,:) = [];
    M2(TF,:) = [];
    
    

    %%
    delta_M = M1 - M2;
    [rho,pval] = corr(M2,M1);
%% loop over upper lower and middle of the u estimation bounds

         NN = size(M1); 
         Number_datapoints = NN(1);
         a = 0.95; % confidence level
         L(1) = 1;% the centre
         L(2) = (Number_datapoints - 1)/chi2inv(((1-a)/2),Number_datapoints-1);% the upper boundary
         L(3) = (Number_datapoints - 1)/chi2inv((1-(1-a)/2),Number_datapoints-1);% the lower boundary
         
         
         diff_L = L(2) - L(1); % this is the percentage of the error
    
        var_M1 = var(M1);
        var_M2 = var(M2);
        var_delta_M = var(delta_M);
    
    
    var_X = 0.5*(var_M1 + var_M2 - var_delta_M);
    var_e1 = 0.5*(var_M1 - var_M2 + var_delta_M);
    var_e2 = 0.5*(var_M2 - var_M1 + var_delta_M);
    
    % centre value
    u_X = (L(1)*var_X)^0.5;
    u_e1 = (L(1)*var_e1)^0.5;
    u_e2 = (L(1)*var_e2)^0.5;
    pu_X = u_X./mean(data.mean_ColumnO3)*100;
    pu_e1 = u_e1./mean( data.mean_ColumnO3)*100;
    pu_e2 = u_e2./mean( data.mean_ColumnO3)*100;
    % error bar
    u_X_err = ((L(2)*var_X)^0.5 - (L(3)*var_X)^0.5)/2;
    u_e1_err = ((L(2)*var_e1)^0.5 - (L(3)*var_e1)^0.5)/2;
    u_e2_err = ((L(2)*var_e2)^0.5 - (L(3)*var_e2)^0.5)/2;
    pu_X_err = u_X_err./mean(data.mean_ColumnO3)*100;
    pu_e1_err = u_e1_err./mean( data.mean_ColumnO3)*100;
    pu_e2_err = u_e2_err./mean( data.mean_ColumnO3)*100;
    
    %% output
    uncertainties = table;
    uncertainties.u_GBS = u_e1;
    uncertainties.u_GBS_err = u_e1_err;
    uncertainties.pu_GBS = pu_e1;
    uncertainties.pu_GBS_err = pu_e1_err;
    uncertainties.u_Brewer = u_e2;
    uncertainties.u_Brewer_err = u_e2_err;
    uncertainties.pu_Brewer = pu_e2;
    uncertainties.pu_Brewer_err = pu_e2_err;
    uncertainties.u_X = u_X;
    uncertainties.u_X_err = u_X_err;
    uncertainties.pu_X = pu_X;
    uncertainties.pu_X_err = pu_X_err;
    uncertainties.rho = rho;
    uncertainties.pval = pval;
    uncertainties.N = NN;
    uncertainties.mean_day = mean(data.day);
    uncertainties.median_day = median(data.day);
    
    
    %linear_fits(M1,M2);
    linear_fits(M2,M1);
    hold all;
    %c = g1.T20(~TF,:);
    %c = data.fd(~TF,:);
    %scatter(M1,M2,20,c,'filled');
    dscatter(M2,M1,'MARKER','o','MSIZE',30,'FILLED',true);
    %scatter(M2,M1,20,c,'filled');
    xlim([-100 100]);
    ylim([-100 100]);
    textbp(['u_G_B_S = ' num2str(u_e1) '; ' num2str(pu_e1) '%']);
    textbp(['u_B_r_e_w_e_r = ' num2str(u_e2) '; ' num2str(pu_e2) '%']);
    textbp(['X = ' num2str(u_X) '; ' num2str(pu_X) '%']);
    N = size(M1);
    textbp(['N = ' num2str(N(1))]);
    xlabel('Brewer ozone [DU]');
    ylabel('GBS ozone [DU]');
    title(fig_title);
    print_setting(1,save_fig,[ fig_title '_scatter_lowf-' num2str(X_Day_step) 'days']);
    %print_setting(1,save_fig,'Brewer_vs_SAOZ_ozone_scatter_filtered');
    
    
    figure;
    h1 = histogram(M1);
    hold on
    h2 = histogram(M2);
    h1.Normalization = 'probability';
    h1.BinWidth = 2;
    h2.Normalization = 'probability';
    h2.BinWidth = 2;
    %xlim([200 550]);
    xlabel('TCO [DU]');
    ylabel('f');
    title(fig_title);
    print_setting(1,save_fig,[ fig_title '_hist_lowf-' num2str(X_Day_step) 'days']);
    %print_setting(1,save_fig,'Brewer_vs_SAOZ_ozone_hist_filtered');
    
end