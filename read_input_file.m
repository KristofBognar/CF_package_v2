function input_table = read_input_file()

input_table = table;

cd ..
cd('CF_package_local');
fid = fopen('CF_input_file.txt','r');
%tline = fgetl(fid);

while ~feof(fid)
    tline = fgetl(fid);
    try
        eval(tline);
    catch
        disp('Warning: one required input might not found! See the following line:')
        disp(tline);
    end
end

