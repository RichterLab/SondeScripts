function [lat_center,lon_center,lat,lon,cach2] = get_track(hurr,sonde_filename)

%From the sonde_filename file, get the day, lat, lon of sonde drop:
fid = fopen(sonde_filename);
C = textscan(fid,'%s','delimiter','\n');
fclose(fid);

%Now pick out relevant data using known character locations:
if (C{1}{16}(7) == '9')
    day = ['19' C{1}{16}(7:12)];
else
    day = ['20' C{1}{16}(7:12)];
end
day = num2str(day);
year = day(1:4);
mon = day(5:6);
day_only = day(7:8);
day = [mon '/' day_only '/' year];

time_char =  C{1}{17}(7:12);
time = str2double(time_char);

cach1=C{1}{17}(7:8); %hour
cach2=str2num(mon)*10000+str2num(day_only)*24+str2num(cach1);

%From the sonde_filename file, get the day, lat, lon of sonde drop:
lat = str2double(C{1}{16}(22:26));
lon = 360 - str2double(C{1}{17}(22:26));

%From the track of the hurricane, get lat and lon as function of time
hurr = hurr(1:end-1); %char
gid = fopen(['./Track_data/', hurr, '.txt']);
D = textscan(gid,'%s','delimiter','\n');
fclose(gid);
% Compare Date, Match Time, and Find Lat/Lon
lat_center=500;
lon_center=600;
for i=4:length(D{1})-1
    %Now pick out relevant data using known character locations:
    date_center = D{1}{i}(1:10);
    time_center = D{1}{i}(12:19);
    hour = time_center(1:2);
    min = time_center(4:5);
    sec = time_center(7:8);
    time_center = [hour min sec];
    time_center = str2double(time_center);
    
    if date_center == day % this should probably be a for loop?
        if time < time_center
            lat_center = str2double(D{1}{i}(25:30));
            lon_center = 360-str2double(D{1}{i}(38:43)); %why?
            break
        end
    end
end
end
