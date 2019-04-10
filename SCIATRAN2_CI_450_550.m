function model = SCIATRAN2_CI_450_550(data_path,simulation_type,plot_path)
%data_path = 'C:\SCIATRAN2_CI\Execute-2.1.5\Data_Out\';
%simulation_type = 'clear_albedo_0dot06_Dec1';
%simulation_type = 'aerosold_lowtran_50km_albedo_0dot06';
%simulation_type = 'cloud_tau12_6to8km_albedo_0dot9';
%simulation_type = 'clear_sky_albedo_0dot9';
fig_size = 0.5; save_fig = 0;
%plot_path = ['H:\work\Eureka\GBS\CI\2010\UTGBS\plots\RTM\' simulation_type];
mkdir(plot_path);
file_nm = 'intensity.dat';
info_nm = 'output_map.inf';
%copyfile(['C:\SCIATRAN2_CI\Execute-2.1.5\Data_Out\' file_nm], plot_path);
%copyfile(['C:\SCIATRAN2_CI\Execute-2.1.5\Data_Out\' info_nm], plot_path);
intensity = importfile([data_path file_nm]);
outputmap = importfile1([data_path info_nm]);

data = intensity;
data(:,1) = [];
data = table2array(data);
data = data';
data = array2table(data);
new_nm = num2str(intensity.wv);
N = size(data);
for i = 1:1:N(2)
    data.Properties.VariableNames{['data' num2str(i)]} = ['wv' new_nm(i,:)];
end

SZA = outputmap.SZA;

%figure; hold all;
% plot(SZA,data.wv360./data.wv550);
% plot(SZA,data.wv405./data.wv550);
plot(SZA,data.wv450./data.wv550);
% plot(SZA,data.wv490./data.wv550);
%legend('360/550 nm','405/550 nm','450/550 nm','490/550 nm');
xlabel('SZA');
ylabel('CI');cd(plot_path);
%print_setting(fig_size, save_fig, ['RTM_CI_vs_SZA' simulation_type '_highlights']);

%% output
model = table;
model.SZA = SZA;
model.CI = data.wv450./data.wv550;

%%
function intensity = importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   INTENSITY = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   INTENSITY = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   intensity = importfile('intensity.dat', 3, 47);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2017/05/25 14:43:42

%% Initialize variables.
if nargin<=2
    startRow = 3;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%13f%23f%23f%23f%23f%23f%23f%23f%23f%23f%23f%23f%23f%23f%23f%23f%23f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
intensity = table(dataArray{1:end-1}, 'VariableNames', {'wv','VarName2','VarName3','VarName4','VarName5','VarName6','VarName7','VarName8','VarName9','VarName10','VarName11','VarName12','VarName13','VarName14','VarName15','VarName16','VarName17','VarName18'});

%%
function outputmap = importfile1(filename, startRow, endRow)
%IMPORTFILE1 Import numeric data from a text file as a matrix.
%   OUTPUTMAP = IMPORTFILE1(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   OUTPUTMAP = IMPORTFILE1(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   outputmap = importfile1('output_map.inf', 5, 21);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2017/05/25 14:57:17

%% Initialize variables.
if nargin<=2
    startRow = 5;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%4f%8f%8f%8f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
outputmap = table(dataArray{1:end-1}, 'VariableNames', {'no','SZA','LOS','AZ','alt'});

