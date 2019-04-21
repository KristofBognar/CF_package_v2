function input_table = read_input_file()

input_table = table;% build empty input information table

cd .. % move to the upper level
cd('CF_package_local'); % move to the local folder, where you have the input file to be read
fid = fopen('CF_input_file.txt','r');% open the local input file

while ~feof(fid) % read evey lines in the file
    tline = fgetl(fid);% get one line of the file
    try
        eval(tline);% okay, the lines in my input file are excutable Matlab format code; for exmaple, "data_path = 'C:/Project/GBS/;'"
    catch
        disp('Warning: one required input might not found! See the following line:')% if it failed to read/excute the line, give a warnning
        disp(tline);% tell me what happend; which line crashed ... 
    end
end

try
    mkdir(input_table.plot_path);% after excute all lines in the "CF_input_file.txt", the input_table should have a column "plot_path", which will be used to save the plots
    status = copyfile('CF_input_file.txt', [input_table.plot_path 'CF_input_file_archive.txt'], 'f');% save the current input settings (the input file) into the plot_path, this is to archive the inputs used for each run/test
    if status == 1 % try to catch some error
        disp('input file has been archived');
    else
        disp('Warning: input file not archived, pls check "read_input_file.m" !');
    end
catch% try to catch some error
    disp('Warning: input file not archived, pls check "read_input_file.m" !');
end