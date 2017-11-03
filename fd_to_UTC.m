function UTC = fd_to_UTC(year,fd)
% this function can convert fractional day of the year to matlab UTC serise
julian = fix(fd);
N = size(fd);
for i = 1:1:N(1)
    [day, month] = Julian2Date(year(i),julian(i));
    hh = (fd(i) - fix(fd(i)))*24;
    mm = (hh - fix(hh))*60;
    HH = fix(hh);
    MM = fix(mm);
    date = [num2str(year(i,:)) '-' num2str(month) '-' num2str(day)];
    time = [num2str(HH) ':' num2str(MM)];
    
    UTC(i,:) = datenum([date ' ' time],'yyyy-mm-dd HH:MM');
end
