function [lat_center,lon_center,lat,lon,time] = get_track(F,sonde_filename)

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

%From the sonde_filename file, get the day, lat, lon of sonde drop:

%THESE ARE THE LAT/LON OF LAUNCH (as opposed to the lat,lon as it falls):
tmpLat = regexp(C{1}{16},'Lat:');
lat = str2double(C{1}{16}(tmpLat+5:tmpLat+11));
tmpLon = regexp(C{1}{17},'Lon:');
lon = str2double(C{1}{17}(tmpLon+5:tmpLon+11));
clear tmpLat tmpLon

lat_center_vec = table2array(F(:,3));
lon_center_vec = table2array(F(:,5));

time_center = datenum(F.Var1(:)+F.Var2(:));

I = find(time<time_center);

if (isempty(I))
    lat_center = NaN;
    lon_center = NaN;
else
   % lat_centers = lat_center_vec(I)
    
    lat_center = lat_center_vec(I(1));
    lon_center = lon_center_vec(I(1));
    
end

end
