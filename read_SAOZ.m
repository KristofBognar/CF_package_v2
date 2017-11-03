function data = read_SAOZ()
% this function read in the SAOZ dSCDs files provided by Andrea, and format
% it to what can be used in CF package
Lat = 80.053;
Lon = -86.416;
Alt = 0.6;
%data_path = 'H:\work\Eureka\SAOZ\SAOZ_2005_2017.mat';
data_path = 'E:\H\work\Eureka\SAOZ\SAOZ_2005_2017_V2.mat';
load(data_path);
C = who('-file',data_path);

N_year = size(C);
%for j = 6:1:N_year(1)
for j = 1:1:11
    eval(['data =' C{j,:} ';']);
    %data = anencours2012;
    %data = data_input;
    try
        data.Year = data.Annee;
    catch
        N_data = size(data);
        vars = C{j,:};
        %data.Annee = repmat(vars(end-3:end),[N_data(1),1]); 
        data.Year = repmat(str2num(vars(end-3:end)),[N_data(1),1]);
    end

data.Fractionalday = data.Jjul;
data.Fractionaltime = (data.Fractionalday - fix(data.Fractionalday))*24;
N = size(data);
eval(['TF' num2str(data.Year(1)) '= false(N(1),1)']);
  
for i = 1:1:N(1)
    %[day, month] = Julian2Date(data.Annee(i),data.Jjul(i));
    [day, month] = Julian2Date(data.Year(i),data.Jjul(i));
    %data.UTC(i) = datenum(data.Annee(i),month,day);
    data.UTC(i) = datenum(data.Year(i),month,day);
    data.DateDDMMYYYY(i)= datetime(datestr(data.UTC(i),'dd-mmm-yyyy'));
    %data.DateDDMMYYYY(i) = datestr(data.UTC(i),'dd-mmm-yyyy');
    [Az El] = SolarAzEl(data.UTC(i),Lat,Lon,Alt);
    data.SolarAzimuthAngle(i) = Az -180;
    
    SZA = 90 - El;
    if abs(SZA - data.SZA(i)) > 0.1
        %disp('warning: check location, calculated SZA do not match SAOZ data');
        %disp(['data SZA :' num2str(data.SZA(i)) ' ; calculated SZA: ' num2str(SZA)]);
        %TF(i) = ture;
        eval(['TF' num2str(data.Year(1)) '(i)= true;']);
%        pause;
    end
    disp_10th = fix(N(1)/10);
    if (i == disp_10th) | (i == 2*disp_10th) | (i == 3*disp_10th) | (i == 4*disp_10th) | (i == 5*disp_10th)| (i == 6*disp_10th)| (i == 7*disp_10th) | (i == 8*disp_10th) | (i == 9*disp_10th)
        disp(['Year' num2str(data.Year(1)) ' : ' num2str(i/N(1)*100) '% finished ...']);
    end
end
disp(['Year ' num2str(data.Year(1)) ' ...']);
data.Elevviewingangle = repmat(90,N(1),1);
data.O3RefZm = zeros(N(1),1);
data.O3ShiftSpectrum = zeros(N(1),1);
data.O3StretchSpectrum1 = zeros(N(1),1);
%data.CI = data.FluxC_450./data.FluxC_550;
data.CI = data.Flux_450./data.Flux_550;
data.O3SlColo3 = data.Q_O3;
data.O3SlErro3 = data.R_O3;
%data.O3RMS = data.ResiduO3Dif;
data.O3RMS = zeros(N(1),1);
data.O4_VIS_293_a203SlColo4 = data.Q_O4_4; 
data.O4_VIS_293_a203SlErro4 = data.R_O4_4;


%figure;
%dscatter(data.SZA,data.CI);

eval(['saoz' num2str(data.Year(1))  '= data;']);
save('SAOZ_formated_temp2.mat');
end;

save('SAOZ_formated_temp2.mat');