function [data,data_output,list_HQ_day] = smooth_flag_v2(data,save_fig)

%load('H:\work\Eureka\GBS\CI\2011\UTGBS\plots\CI_callibration\u1_2011_cal_CI.mat');
plot_residual_statistic = 0;

%data = data_output; % note will use un-filtered data as read in
data_output = data;
%% filter
TF_SZA = data.SZA <=92;
data(~TF_SZA,:) = [];

%% smooth CI
data.UTC = datenum(data.DateDDMMYYYY) + data.Fractionaltime./24;
data.LT = data.UTC - 5/24;
data.LT_DateDDMMYYY = datestr(data.LT,'dd-mmm-yyyy');
data.LT_DateDDMMYYY = datetime(data.LT_DateDDMMYYY,'InputFormat','dd-MMM-yyyy');
N_days = unique(data.LT_DateDDMMYYY);
N = size(N_days);
N_total = size(data);
figure;hold all;
data.TF_smooth = zeros(N_total(1),1);
data.TF_smooth_Gle = zeros(N_total(1),1);
data.TF_smooth_Wag = zeros(N_total(1),1);
for i =1:1:N(1)

    TF = data.LT_DateDDMMYYY == N_days(i);
    x = data.Fractionalday(TF);
    y = data.cal_CI(TF);
    plot(x,y,'.');
    try
    yy = smooth(x,y,0.5,'rloess');
    plot(x,yy);
    TF_smooth = abs((y - yy)./yy) > 0.1;
    catch
        TF_smooth = false;
    end
    data.TF_smooth(TF) = TF_smooth;
    
   %% Glelen's method
   %x = data.Fractionalday(TF);
   %x = data.SZA(TF);
   x = data.Fractionaltime(TF);
   y = data.cal_CI(TF);
   %beta0 = [0,0,0,0,0,0,0];
   beta0 = [1.5,0.1,1,10,0.05,3,10];
   modelfun = @(b,x)(b(1) + b(2)*sin(b(3)*x-b(4)) + b(5)*sin(b(6)*x-b(7)));
   opts = statset('nlinfit');
   opts.RobustWgtFun = 'bisquare';
   beta1 = nlinfit(x,y,modelfun,beta0,opts);
   beta = nlinfit(x,y,modelfun,beta1,opts);
   yy_Gle = (beta(1) + beta(2)*sin(beta(3)*x-beta(4)) + beta(5)*sin(beta(6)*x-beta(7)));  
   plot(data.Fractionalday(TF),yy_Gle);
   TF_smooth_Gle = abs((y - yy_Gle)./yy_Gle) > 0.1;
   data.TF_smooth_Gle(TF) = TF_smooth_Gle;
   
   
   %% Wagner's method
   N_meas = size(data.Fractionaltime(TF));
   x = data.Fractionaltime(TF);
   y = data.cal_CI(TF);
   TF_smooth_Wag = zeros(N_meas(1),1) ;
   for j = 2:1:N_meas(1)-1
    t1 = x(j) - x(j-1);
    t2 = x(j+1) - x(j);
    y_n0 = y(j-1);
    y_n1 = y(j);
    y_n2 = y(j+1);
    if ~(t1 > 0.5 | t2 > 0.5)
        TF_smooth_Wag(j,1) = 2*((t1*y_n2 + t2*y_n0)/(t1*t2*(t1+t2)) - y_n1/(t1*t2));
    end
   end
   data.TF_smooth_Wag(TF) = TF_smooth_Wag;
   
   
   disp([num2str(i/N(1)*100) ' % finished ...']);
end

data.TF_smooth = logical(data.TF_smooth);
data.TF_smooth_Gle = logical(data.TF_smooth_Gle);
TF = abs(data.TF_smooth_Wag) > 10;
data.TF_smooth_Wag_1 = TF;


plot(data.Fractionalday(data.TF_smooth),data.cal_CI(data.TF_smooth),'rx');
plot(data.Fractionalday(data.TF_smooth_Gle),data.cal_CI(data.TF_smooth_Gle),'bx');
plot(data.Fractionalday(data.TF_smooth_Wag_1),data.cal_CI(data.TF_smooth_Wag_1),'gx');

p = sum(data.TF_smooth)./N_total(1)*100;
p_Gle = sum(data.TF_smooth_Gle)./N_total(1)*100;
p_Wag_1 = sum(data.TF_smooth_Wag_1)./N_total(1)*100;
disp([num2str(p) ' % data marked as high TSI']);
disp([num2str(p_Gle) ' % data marked as high TSI (Gielen method)']);
disp([num2str(p_Wag_1) ' % data marked as high TSI (Wagner method)']);
textbp([num2str(p) ' % data marked as high TSI (red x)']);
textbp([num2str(p_Gle) ' % data marked as high TSI (Gielen method, blue x)']);
textbp([num2str(p_Wag_1) ' % data marked as high TSI (Wagner method, green x)']);
xlabel('Fractional day of year');
ylabel('Calibrated CI');
print_setting(1,save_fig,['CI smoothness1']);

figure; hold all;
gscatter(data.SZA,data.cal_CI,data.clear);
plot(data.SZA(data.TF_smooth),data.cal_CI(data.TF_smooth),'rx');
plot(data.SZA(data.TF_smooth_Gle),data.cal_CI(data.TF_smooth_Gle),'bx');
plot(data.SZA(data.TF_smooth_Wag_1),data.cal_CI(data.TF_smooth_Wag_1),'gx');
xlabel('SZA');
ylabel('Calibrated CI');
print_setting(1,save_fig,['CI smoothness2']);


%% smooth O4
% data.UTC = datenum(data.DateDDMMYYYY) + data.Fractionaltime./24;
% data.LT = data.UTC - 5/24;
% data.LT_DateDDMMYYY = datestr(data.LT,'dd-mmm-yyyy');
% data.LT_DateDDMMYYY = datetime(data.LT_DateDDMMYYY,'InputFormat','dd-MMM-yyyy');
% N_days = unique(data.LT_DateDDMMYYY);
% N = size(N_days);
% N_total = size(data);
figure;hold all;
data.TF_smooth_O4 = zeros(N_total(1),1);
data.TF_smooth_O4_Gle = zeros(N_total(1),1);
data.TF_smooth_O4_Wag = zeros(N_total(1),1);
O4VCD_offset = 1.41e43;
for i =1:1:N(1)

    TF = data.LT_DateDDMMYYY == N_days(i);
    x = data.Fractionalday(TF);
    %y = data.O4_VIS_293_a203SlColo4(TF);
    %y = data.O4_VIS_293_a203SlColo4(TF) + O4VCD_offset;
    %y_err = data.O4_VIS_293_a203SlErro4(TF);
    try
        y = data.O3SlColo4(TF) + O4VCD_offset;
        y_err = data.O3SlErro4(TF);
    catch
        y = data.Q_O4_4(TF) + O4VCD_offset;
        y_err = data.R_O4_4(TF);
    end
    
    plot(x,y,'.');
    try
        yy = smooth(x,y,0.5,'rloess');
        plot(x,yy);
        TF_smooth_O4 = abs((y - yy)./yy) > 0.2;
    catch
        TF_smooth_O4 = false;
    end
    upper_err = y + y_err;
    lower_err = y - y_err;
    TF_err = (upper_err >= 0) & (lower_err <= 0);
    TF_smooth_O4 = TF_smooth_O4 & ~(TF_err);
    data.TF_smooth_O4(TF) = TF_smooth_O4;
    print_setting(1,save_fig,['O4 smoothness1']);
    
%    %% Glelen's method
%    x = data.Fractionaltime(TF);
%    y = data.O4_VIS_293_a203SlColo4(TF);
%    beta0 = [5e41,0.1e43,1,10,0.05e43,3,10];
%    modelfun = @(b,x)(b(1) + b(2)*sin(b(3)*x-b(4)) + b(5)*sin(b(6)*x-b(7)));
%    opts = statset('nlinfit');
%    opts.RobustWgtFun = 'bisquare';
%    beta1 = nlinfit(x,y,modelfun,beta0,opts);
%    beta = nlinfit(x,y,modelfun,beta1,opts);
%    yy_Gle = (beta(1) + beta(2)*sin(beta(3)*x-beta(4)) + beta(5)*sin(beta(6)*x-beta(7)));  
%    plot(data.Fractionalday(TF),yy_Gle);
%    TF_smooth_O4_Gle = abs((y - yy_Gle)./yy_Gle) > 0.2;
%    data.TF_smooth_O4_Gle(TF) = TF_smooth_O4_Gle;
   
   
   %% Wagner's method
   N_meas = size(data.Fractionaltime(TF));
   x = data.Fractionaltime(TF);
   %y = data.O4_VIS_293_a203SlColo4(TF);
   TF_smooth_O4_Wag = zeros(N_meas(1),1) ;
   for j = 2:1:N_meas(1)-1
    t1 = x(j) - x(j-1);
    t2 = x(j+1) - x(j);
    y_n0 = y(j-1);
    y_n1 = y(j);
    y_n2 = y(j+1);
    if ~(t1 > 0.5 | t2 > 0.5)
        TF_smooth_O4_Wag(j,1) = 2*((t1*y_n2 + t2*y_n0)/(t1*t2*(t1+t2)) - y_n1/(t1*t2));
    end
   end
   data.TF_smooth_O4_Wag(TF) = TF_smooth_O4_Wag;
   
   disp([num2str(i/N(1)*100) ' % O4 smoothness finished ...']);
   %disp([num2str(i/N(1)*100) ' % finished ...']);
end

data.TF_smooth_O4 = logical(data.TF_smooth_O4);
% data.TF_smooth_O4_Gle = logical(data.TF_smooth_O4_Gle);
TF = abs(data.TF_smooth_O4_Wag) > 2e44;
data.TF_smooth_O4_Wag_1 = TF;


%% output data (note, data_output retained all measurements (all SZA)!)
N_raw = size(data_output);


%data.HQ_index = data.clear & ~data.TF_smooth & ~data.TF_smooth_O4;
%data.HQ_index = ~data.clear | data.TF_smooth | data.TF_smooth_O4;
data.HQ_index = data.cloudy | data.TF_smooth | data.TF_smooth_O4;

data_output.TF_smooth_alter = logical(zeros(N_raw(1),1));
data_output.TF_smooth_alter_O4 = logical(zeros(N_raw(1),1));
data_output.HQ_index_alter = logical(zeros(N_raw(1),1));
data_output.TF_smooth_alter(TF_SZA,:) = data.TF_smooth;
data_output.TF_smooth_alter_O4(TF_SZA,:) = data.TF_smooth_O4;
data_output.HQ_index_alter(TF_SZA,:) = data.HQ_index;



%% plots

%y = data.O4_VIS_293_a203SlColo4 + O4VCD_offset;
try
    y = data.O3SlColo4 + O4VCD_offset;
catch
    y = data.Q_O4_4 + O4VCD_offset;
end
plot(data.Fractionalday(data.TF_smooth_O4),y(data.TF_smooth_O4),'rs');
% plot(data.Fractionalday(data.TF_smooth_O4_Gle),data.O4_VIS_293_a203SlColo4(data.TF_smooth_O4_Gle),'bs');
plot(data.Fractionalday(data.TF_smooth_O4_Wag_1),y(data.TF_smooth_O4_Wag_1),'gs');

p = sum(data.TF_smooth_O4)./N_total(1)*100;
% p_Gle = sum(data.TF_smooth_O4_Gle)./N_total(1)*100;
p_Wag_1 = sum(data.TF_smooth_O4_Wag_1)./N_total(1)*100;
disp([num2str(p) ' % data marked as high O_4 TSI']);
% disp([num2str(p_Gle) ' % data marked as high O_4 TSI (Gielen method)']);
disp([num2str(p_Wag_1) ' % data marked as high O_4 TSI (Wagner method)']);
textbp([num2str(p) ' % data marked as high O_4 TSI (red s)']);
% textbp([num2str(p_Gle) ' % data marked as high O_4 TSI (Gielen method, blue s)']);
textbp([num2str(p_Wag_1) ' % data marked as high O_4 TSI (Wagner method, green s)']);
xlabel('Fractional day of year');
ylabel('O_4 dSCD');
print_setting(1,save_fig,['O4 smoothness1']);


figure; hold all;
gscatter(data.SZA,data.cal_CI,data.clear);
plot(data.SZA(data.TF_smooth),data.cal_CI(data.TF_smooth),'rx');
% plot(data.SZA(data.TF_smooth_Gle),data.cal_CI(data.TF_smooth_Gle),'bx');
plot(data.SZA(data.TF_smooth_Wag_1),data.cal_CI(data.TF_smooth_Wag_1),'gx');
plot(data.SZA(data.TF_smooth_O4),data.cal_CI(data.TF_smooth_O4),'rs');
% plot(data.SZA(data.TF_smooth_O4_Gle),data.cal_CI(data.TF_smooth_O4_Gle),'bs');
plot(data.SZA(data.TF_smooth_O4_Wag_1),data.cal_CI(data.TF_smooth_O4_Wag_1),'gs');
xlabel('SZA');
ylabel('Calibrated CI');
print_setting(1,save_fig,['CI marked']);

% %% smooth O4
% N_days = unique(data.DateDDMMYYYY);
% N = size(N_days);
% N_total = size(data);
% figure;hold all;
% data.TF_O4_smooth= zeros(N_total(1),1);
% for i =1:1:N(1)
%    TF = data.DateDDMMYYYY == N_days(i);
%    plot(data.Fractionalday(TF),data.O4_VIS_293_a203SlColo4(TF),'.');
%    yy = smooth(data.Fractionalday(TF),data.O4_VIS_293_a203SlColo4(TF),0.5,'rloess');
%    plot(data.Fractionalday(TF),yy);
%    TF_O4_smooth = abs((data.O4_VIS_293_a203SlColo4(TF) - yy)./yy) > 0.2;
%    data.TF_O4_smooth(TF) = TF_O4_smooth;
% end
% data.TF_O4_smooth = logical(data.TF_O4_smooth);
% plot(data.Fractionalday(data.TF_O4_smooth),data.O4_VIS_293_a203SlColo4(data.TF_O4_smooth),'rx');
% p = sum(data.TF_O4_smooth)./N_total(1)*100;
% disp([num2str(p) ' % data marked as high O4 TSI']);
% textbp([num2str(p) ' % data marked as high O_4 TSI (red x)']);
% xlabel('Fractional day of year');
% ylabel('O_4 dSCDs');

figure; hold all;
gscatter(data.SZA,data.cal_CI,data.clear);
plot(data.SZA(data.TF_smooth),data.cal_CI(data.TF_smooth),'rx');
plot(data.SZA(data.TF_smooth_O4),data.cal_CI(data.TF_smooth_O4),'bx');
xlabel('SZA');
ylabel('Calibrated CI');
print_setting(1,save_fig,['CI grouped month']);

figure;hold all;
plot(data.Fractionalday(~data.HQ_index),data.cal_CI(~data.HQ_index),'.');
plot(data.Fractionalday(data.HQ_index),data.cal_CI(data.HQ_index),'.');
plot(data.Fractionalday(data.TF_smooth),data.cal_CI(data.TF_smooth),'s');
plot(data.Fractionalday(data.TF_smooth_O4),data.cal_CI(data.TF_smooth_O4),'s');
legend('high quality data', 'removed data','filtered by CI TSI','filtered by O_4 TSI');
p_hq = (1 - sum(data.HQ_index)/N_total(1))*100;
textbp([num2str(p_hq) ' % data flaged as HQ']);
xlabel('Day of the year');
ylabel('Calibrated CI');
print_setting(1,save_fig,['CI marks 2']);


for j = 1:1:N(1)
    TF = data.LT_DateDDMMYYY == N_days(j);
    TF_HQ_day(j) = sum(data.HQ_index(TF)) == 0;
end
list_HQ_day = N_days(TF_HQ_day);
p_hq_day = sum(TF_HQ_day)/N(1)*100;
textbp([num2str(p_hq_day) '% (' num2str(sum(TF_HQ_day)) ') days (out of' num2str(N(1)) ') flaged as HQ']);

%% residual statistic
if plot_residual_statistic == 1;
    
    bin = 2e-5;
    TF = data.clear == 1;
    figure;hold all;
    h1=histogram(data.O3293KRMS(TF));
    h1.BinWidth = bin;
    h1.Normalization = 'probability';
    TF = data.mediocre == 1;
    h3=histogram(data.O3293KRMS(TF));
    h3.BinWidth = bin;
    h3.Normalization = 'probability';
    TF = data.cloudy == 1;
    h2=histogram(data.O3293KRMS(TF));
    h2.BinWidth = bin;
    h2.Normalization = 'probability';
    legend('clear','mediocre','cloudy');
    xlim([0 1.5e-3]);
    ylabel('f');
    xlabel('O3 RMS');
    
    %TF = data.O3293Kprocessing_error == 1;
    TF = data.NO2298KRMS > 1e-2;
    data(TF,:) = [];
    bin = 2e-5;
    TF = data.clear == 1;
    figure;hold all;
    h1=histogram(data.NO2298KRMS(TF));
    h1.BinWidth = bin;
    h1.Normalization = 'probability';
    TF = data.mediocre == 1;
    h3=histogram(data.NO2298KRMS(TF));
    h3.BinWidth = bin;
    h3.Normalization = 'probability';
    TF = data.cloudy == 1;
    h2=histogram(data.NO2298KRMS(TF));
    h2.BinWidth = bin;
    h2.Normalization = 'probability';
    legend('clear','mediocre','cloudy');
    xlim([0 2e-3]);
    ylabel('f');
    xlabel('NO2 RMS');
    
    
    
    
    bin = 2e-5;
    TF = data.clear == 1;
    figure;hold all;
    h1=histogram(data.O4_VIS_293_a203RMS(TF));
    h1.BinWidth = bin;
    h1.Normalization = 'probability';
    TF = data.mediocre == 1;
    h3=histogram(data.O4_VIS_293_a203RMS(TF));
    h3.BinWidth = bin;
    h3.Normalization = 'probability';
    TF = data.cloudy == 1;
    h2=histogram(data.O4_VIS_293_a203RMS(TF));
    h2.BinWidth = bin;
    h2.Normalization = 'probability';
    legend('clear','mediocre','cloudy');
    xlim([0 2e-3]);
    ylabel('f');
    xlabel('O4 RMS');
end



