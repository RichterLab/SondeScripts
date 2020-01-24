function [lat_center,lon_center,lat,lon,cach2] = get_track(hurr,sonde_filename)

%From the sonde_filename file, get the day, lat, lon of sonde drop:
fid = fopen(sonde_filename);
C = textscan(fid,'%s','delimiter','\n');
fclose(fid);

%Now pick out relevant data using known character locations:
if (C{1}{16}(7) == '9')
    date_sonde = ['19' C{1}{16}(7:12)];
else
    date_sonde = ['20' C{1}{16}(7:12)];
end
year = date_sonde(1:4);
mon = date_sonde(5:6);
day = date_sonde(7:8);

time_sonde =  C{1}{17}(7:12);

hr = time_sonde(1:2);
min = time_sonde(3:4);
sec = time_sonde(5:6);

%Get a numeric value of date/time:
time = datenum(str2double(year),str2double(mon),str2double(day),str2double(hr),str2double(min),str2double(sec));

%Format sonde time output in same way as done for EBT data (i.e. RMW)
cach1=C{1}{17}(7:8); %hour
cach2=str2double(mon)*10000+str2double(day)*24+str2double(cach1);

%From the sonde_filename file, get the day, lat, lon of sonde drop:

%THESE ARE THE LAT/LON OF LAUNCH (as opposed to the lat,lon as it falls):
lat = str2double(C{1}{16}(22:26));
lon = 360 - str2double(C{1}{17}(22:26));


%From the track of the hurricane, get lat and lon as function of time
hurr = hurr(1:end-1); %char
gid = fopen(['./Track_data/', hurr, '.txt']);
D = textscan(gid,'%s','delimiter','\n');
fclose(gid);

% Compare Date, Match Time, and Find Lat/Lon

%To hopefully catch (via breaking) in case a center comparison isn't found
lat_center=NaN;
lon_center=NaN;

for i=4:length(D{1})-1
    %Now pick out relevant data using known character locations:
    date_center = D{1}{i}(1:10);
    time_center = D{1}{i}(12:19);
    year_center = date_center(7:10);
    mon_center = date_center(1:2);
    day_center = date_center(4:5);
    hour = time_center(1:2);
    min = time_center(4:5);
    sec = time_center(7:8);
    time_center = datenum(str2double(year_center),str2double(mon_center),str2double(day_center),str2double(hour),str2double(min),str2double(sec));
    
    if time < time_center
        lat_center = str2double(D{1}{i}(25:30));
        lon_center = 360-str2double(D{1}{i}(38:43));
        break
    end
end

end
