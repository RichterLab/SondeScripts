%%
clear all
close all
clc

%Split into 10m/s-wide bins from 0 to 80 m/s (each with its own mean profile)
%Split vertical into 5m-wide bins

cpa = 1.005;  %kJ/kg-K
cpl = 4.19;   %kJ/kg-K
cpv = 1.86;   %kJ/kg-K
Lv = 2260;    %kJ/kg
rho = 1.2;   %kg/m^3

% Split data into bins for height, temp., windspeed, and relative humidity.
% Focus only on breaking up height bins first.

pbl_height = 10000;  %Height over which mean is taken (Powell uses 500m)
max_height = 10000;  %Height over which the fit is done
min_height = 10;   %Bottom height over which fit is done
height_interval = 10;
num_z = (max_height-min_height)/height_interval;

% Limits of wind speed bins
max_wind = 100;
wind_interval = 10;
num_wind_bins = max_wind/wind_interval;

%Limits of radius bins -- in normalized units (R/RMW)
max_rad = 12;
rad_interval = 1;
num_rad_bins = max_rad/rad_interval;

%According to Guiquan, do it in more discrete intervals:
%num_rad_bins = 4; % in the eyewall or outside the eyewall
%0-0.5; 0.1-1; 1-2; 2-3; 3-4; >4
%0-1; 1-2; 2-3; >3

% Making blank matrices to save spaces for data
% [height_bin, radii_bin, wind_bin]
mean_U_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
mean_Ur_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
mean_Ut_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
numvecU = zeros(num_z,num_rad_bins,num_wind_bins);
numsonde = zeros(num_rad_bins,num_wind_bins);
mean_T_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
% numvecT = zeros(num_z,num_rad_bins,num_theta_bins);
mean_RH_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
mean_zU_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
mean_z_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
mean_p_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
mean_q_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
mean_k_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
mean_theta_profiles = zeros(num_z,num_rad_bins,num_wind_bins);

all_U_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_T_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_RH_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_zU_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_z_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_p_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_q_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_k_profiles = cell(num_z,num_rad_bins,num_wind_bins);
all_theta_profiles = cell(num_z,num_rad_bins,num_wind_bins);
%Guiquan ssttosave = cell(num_rad_bins,num_theta_bins);

zplot = min_height+0.5*height_interval:height_interval:max_height-0.5*height_interval;
rplot = 0.5*rad_interval:rad_interval:max_rad-0.5*rad_interval;
Uplot = 0.5*wind_interval:wind_interval:max_wind-0.5*wind_interval;

%Without basically empty ones (as determined by composite wind height-radius plots):
% hurrvec = { ...
%     'Danny1997/' 'Erika1997/' 'Guillermo1997/' ...
%     'Bonnie1998/' 'Danielle1998/' 'Earl1998/' 'Georges1998/' 'Mitch1998/' ...
%     'Bret1999/' 'Dennis1999/' 'Dora1999/' 'Emily1999/' 'Eugene1999/' 'Floyd1999/' 'Irene1999/' ...
%     'Helene2006/' ...
%     'Felix2007/' 'Ingrid2007/' 'Karen2007/' 'Noel2007/' ...
%     'Dolly2008/' 'Fay2008/' 'Gustav2008/' 'Hanna2008/' 'Ike2008/' 'Kyle2008/' 'Paloma2008/' ...
%     'Ana2009/' 'Bill2009/' 'Danny2009/' ...
%     'Alex2010/' 'Earl2010/' 'Karl2010/' 'Richard2010/' 'Tomas2010/' ...
%     'Dora2011/' 'Don2011/' 'Emily2011/' 'Irene2011/' 'Katia2011/' 'Rina2011/' ...
%     'Isaac2012/' 'Leslie2012/' 'Sandy2012/' ...
%     'Gabrielle2013/' 'Ingrid2013/' 'Karen2013/' ...
%     'Arthur2014/' 'Bertha2014/' 'Cristobal2014/' 'Edouard2014/' 'Gonzalo2014/' 'Iselle2014/' 'Julio2014/' 'Simon2014/' 'Ana2014/' ...
%     'Bill2015/' 'Carlos2015/' 'Danny2015/' 'Ela2015/' 'Erika2015/' 'Guillermo2015/' 'Hilda2015/' 'Joaquin2015/' 'Kate2015/' 'Oho2015/' 'Patricia2015/' ...
%     'Bonnie2016/' 'Colin2016/' 'Danielle2016/' 'Darby2016/' 'Earl2016/' 'Hermine2016/' 'Javier2016/' 'Julia2016/' 'Karl2016/' 'Lester2016/' 'Matthew2016/' 'Nicole2016/' 'Otto2016/' ...
%     'Franklin2017/' 'Harvey2017/' 'Irma2017/' 'Jose2017/' 'Katia2017/' 'Maria2017/' 'Nate2017/' 'Norma2017/'};

% Without ones with no profiles used in the processing
 hurrvec = { ... 
      'Danny1997/' 'Erika1997/' ...
      'Bonnie1998/' 'Danielle1998/' 'Earl1998/' 'Georges1998/' 'Mitch1998/' ...
      'Bret1999/' 'Dennis1999/' 'Floyd1999/' 'Irene1999/' ...
      'Edouard2002/' 'Fay2002/' 'Gustav2002/' 'Hanna2002/' 'Isidore2002/' 'Lili2002/' ...
      'Claudette2003/' 'Erika2003/' 'Fabian2003/' 'Isabel2003/' 'Odette2003/' ...
      'Bonnie2004/' 'Charley2004/' 'Frances2004/' 'Ivan2004/' 'Jeanne2004/' ...
      'Arlene2005/' 'Cindy2005/' 'Dennis2005/' 'Emily2005/' 'Franklin2005/' 'Gert2005/' 'Irene2005/' 'Katrina2005/' 'Ophelia2005/' 'Rita2005/' 'Wilma2005/' 'Gamma2005/' ...
      'Ernesto2006/' 'Helene2006/' ...
      'Felix2007/' 'Ingrid2007/' 'Karen2007/' ...
      'Dolly2008/' 'Fay2008/' 'Gustav2008/' 'Ike2008/' 'Kyle2008/' 'Paloma2008/' ...
      'Bill2009/' ...
      'Alex2010/' 'Earl2010/' 'Karl2010/' 'Richard2010/' 'Tomas2010/' ...
      'Dora 2011/' 'Don2011/' 'Emily2011/' 'Irene2011/' 'Katia2011/' 'Rina2011/' ...
      'Isaac2012/' 'Leslie2012/' 'Sandy2012/' ...
      'Gabrielle2013/' 'Ingrid2013/' 'Karen2013/' ...
      'Arthur2014/' 'Bertha2014/' 'Cristobal2014/' 'Edouard2014/' 'Gonzalo2014/' 'Iselle2014/' 'Julio2014/' 'Simon2014/' ...
      'Carlos2015/' 'Danny2015/' 'Erika2015/' 'Guillermo2015/' 'Hilda2015/' 'Joaquin2015/' 'Kate2015/' 'Patricia2015/' ...
      'Bonnie2016/' 'Danielle2016/' 'Earl2016/' 'Hermine2016/' 'Javier2016/' 'Julia2016/' 'Karl2016/' 'Matthew2016/' 'Nicole2016/' 'Otto2016/' ...
      'Franklin2017/' 'Harvey2017/' 'Irma2017/' 'Jose2017/' 'Katia2017/' 'Maria2017/' 'Nate2017/' 'Norma2017/' ...
      'Alberto2018/' 'Chris2018/' 'Florence2018/' 'Gordon2018/' 'Hector2018/' 'Isaac2018/' 'Michael2018/' 'Olivia2018/'
};

%%% Load both TCOBS and EBT RMW data, use TCOBS if available, EBT if not

%%% 1. RMW from Dan Chavas, every 6 hours%%%%%
% load Atlantic and Pacific data separately then combine 
atl = load('./RWS/ebtrk_atlc_1988_2018.mat','Year_EBT','Month_EBT','Day_EBT','HourUTC_EBT','StormName_EBT','rmkm_EBT','Vmms_EBT','LatN_EBT','LonW_EBT');
pac = load('./RWS/ebtrk_epac_2001_2018.mat','Year_EBT','Month_EBT','Day_EBT','HourUTC_EBT','StormName_EBT','rmkm_EBT','Vmms_EBT','LatN_EBT','LonW_EBT');

Year_EBT = [atl.Year_EBT; pac.Year_EBT];
Month_EBT = [atl.Month_EBT; pac.Month_EBT];
Day_EBT = [atl.Day_EBT; pac.Day_EBT];
HourUTC_EBT = [atl.HourUTC_EBT; pac.HourUTC_EBT];
StormName_EBT = [atl.StormName_EBT; pac.StormName_EBT];
rmkm_EBT = [atl.rmkm_EBT; pac.rmkm_EBT];
Vmms_EBT = [atl.Vmms_EBT; pac.Vmms_EBT];
LatN_EBT = [atl.LatN_EBT; pac.LatN_EBT];
LonW_EBT = [atl.LonW_EBT; pac.LonW_EBT];
rmw_name_EBT = strcat(StormName_EBT,num2str(Year_EBT),'/');
rmw_time_EBT = datenum([Year_EBT';Month_EBT';Day_EBT';HourUTC_EBT';zeros(size(Day_EBT))';zeros(size(Day_EBT))']'); % time in matlab datenum format
clear atl pac

%%% 2. RMW from Jonathan Vigh, every 1 hour%%%%%
%%% Attention: The datatype is different than 'ebtrk_atlc_1988_2017.mat'%%%
load('./RWS/TCOBS_0.42_20180209.nc.mat');

cach = num2str(ATCF_stormyear); %extract year added to rmw_name
rmw_name_TC = strcat(ATCF_stormname,cach,'/');
convDist = 1.852;   % Conversion factor for naut miles to km
convSpeed = 0.51444; % COnversion factor from knots to m/s
rmw_time_TC = nan(size(TCOBS_RMW));
for currStorm =1:length(rmw_name_TC)      % Cycle through storms
    timeStart = BT_yyyymmddhhmm{1,currStorm};
    tcTstart = datenum([str2num(timeStart(1:4)) str2num(timeStart(5:6)) str2num(timeStart(7:8)) ...
        str2num(timeStart(9:10)) str2num(timeStart(11:12)) 00]);
    numObs = sum(~isnan(TCOBS_time_since_first_data_point(:,currStorm)));
    for w=1:numObs   % Cycle through TCOBS for that storm
        rmw_time_TC(w,currStorm) = datenum([str2num(timeStart(1:4)) str2num(timeStart(5:6)) str2num(timeStart(7:8)) ...
            str2num(timeStart(9:10))+TCOBS_time_since_first_data_point(w,currStorm) str2num(timeStart(11:12)) 00]);
    end 
end
rmkm_TC = convDist.*TCOBS_RMW;
Vmms_TC = convSpeed.*TCOBS_Vmax;

clearvars cach

count_hurr=0;
count_sonde = 1;
for hurr = hurrvec % 1. loop every hurricane
    count_hurr=count_hurr+1;
    mean_U_profiles_onestorm = zeros(num_z,num_rad_bins,num_wind_bins);
    numvecU_onestorm = zeros(num_z,num_rad_bins,num_wind_bins);
    
    fprintf('Starting storm: %s\n',hurr{1});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check if storm is in TCOBS or EBT
    clear isTCOBS isEBT
    isTCOBS = ismember(upper(hurr),rmw_name_TC);
    if isTCOBS<1
       isEBT = ismember(upper(hurr),rmw_name_EBT);
    else
        isEBT = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (isEBT + isTCOBS > 0) %(~isnan(cach)) % ensure hurr name of dropsondes corresponds to a rmw_name from the track data
        if isTCOBS>0
            clearvars rmw_time rmw_rad_hurr
            cach = strmatch(upper(hurr), rmw_name_TC); % RMW location of specific hurricane
            cach1 = rmw_time_TC(:,cach); % RMW time of specific hurricane
            cach2 = rmkm_TC(:,cach); % RMW radius of specific hurricane
            cach3 = Vmms_TC(:,cach); % RMW wind speed of specific hurricane
            rmw_time_hurr = cach1(cach1>0&cach2>0);%remove nan
            rmw_rad_hurr  = cach2(cach1>0&cach2>0);%remove nan
            rmw_speed_hurr= cach3(cach1>0&cach2>0);%remove nan
        elseif isEBT>0
            clearvars rmw_time rmw_rad_hurr
            cach = strmatch(upper(hurr), rmw_name_EBT);
            cach1 = rmw_time_EBT(cach); % RMW time of specific hurricane
            cach2 = rmkm_EBT(cach); % RMW radius of specific hurricane
            cach3 = Vmms_EBT(cach); % RMW wind speed of specific hurricane
            rmw_time_hurr = cach1(cach1>0&cach2>0);%remove nan
            rmw_rad_hurr  = cach2(cach1>0&cach2>0);%remove nan
            rmw_speed_hurr= cach3(cach1>0&cach2>0);%remove nan
        end
       
        clearvars cach1 cach2 cach3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Before cycling through the sondes, open and read this storm's track file:
        trackfile = strcat('./Track_data/',hurr{1}(1:end-1),'.txt');
        F = readtable(trackfile,'Format','%{MM/dd/yyyy}D %{hh:mm:ss}T %f %s %f %s');
        
        filedirtmp = strcat('./all_storms/',hurr);
        filedir = filedirtmp{1};
        files = dir([filedir '*.frd']);
        
        %Try loading all frd files from a single storm, stored using 'read_frd_files.m' to reduce reading time of each individual sonde
        frd_filename = strcat('./all_storms/',hurr{1}(1:end-1),'_frd_files.mat');
        load(frd_filename);
        
        for i=1:length(files) % 2. loop every dropsonde
            
            %fprintf('Opening file: %s\n',[filedir files(i).name]);
            %datstruct = importdata([filedir files(i).name],' ',21);
            %dat = datstruct.data;
            dat = frddat{i};
            
            %Assign the values - only take data where z,temp,p,WS,RH above 0
            ztmp = dat(:,6); % [m]
            Ttmp = dat(:,4); % [C]
            ptmp = dat(:,3); % [mb]
            WStmp = dat(:,8); % [m/s]
            RHtmp = dat(:,5); % [%]
            Utmp = dat(:,9); % [m/s] -- the zonal velocity component
            Vtmp = dat(:,10); % [m/s] -- the meridonal velocity component
            lattmp = dat(:,18); %Deg N
            lontmp = dat(:,19); %Deg E
            
            %Clean the data - if any of them have a negative signal, remove for
            %all profiles (only want a profile if all data is there)
            % This is the set I'll be plotting for now (all with full at end)
            z = ztmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0);
            T = Ttmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0) + 273;
            p = ptmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0);
            WS = WStmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0) ;
            U = Utmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0);
            V = Vtmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0);
            RH = RHtmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0);
            latvec = lattmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0);
            lonvec = lontmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0);
            Ut = zeros(size(U));
            Ur = zeros(size(U));
            
            
            if (~isempty(z))  %Only do the computations if data is left after the clean
                
                %Compute potential temperature:
                psurf = p(length(p));
                theta = T.*(psurf./p).^(2/7);
                
                %Compute the specific humidity q:
                
                %Need saturation vapor pressure:
                %Formula from Dingman's hydrology book, eq. D-7 ("Magnus relation")
                %Needs T in celsius, gives estar in mb (hPa)
                estar = 6.11*exp(17.3*(T-273)./((T-273)+237.3));
                e = (RH/100).*estar;
                
                %p is in mb as well, so can divide:
                eoverp = e./p;
                q = 0.622*eoverp./(1-0.378*eoverp);
                
                %Now compute enthalpy k, as defined in Emanuel (1995):
                %k will have units of kJ/kg
                
                %Test using temperature:
                k = (cpa*(1-q) + cpl*q).*theta + Lv*q;
                
                %Compute mean velocity below pbl_height (Powell uses 500m):
                %WSmean = mean(WS(z<pbl_height));
                
                % Try using the maximum wind speed below 100 m
                WSmean = nanmax(WS(z<100));
                
                %Trying by using the lowest recorded wind speed rather than a PBL average
%                 if (z(end) < 100)  %Only do this if the lowest recording is beneath 100m
%                     WSmean = WS(end);
%                 else
%                     WSmean = NaN;
%                 end
                
%                 %Testing old way of binning by wind speed:
%                 WSmean = mean(WS(z<50));
%                 
%                 %Trying to bin by the wind speed near 10m
                idx_test = find(z<10,1);
                if (idx_test > 0)
                    WSmean = WS(idx_test);
                else
                    WSmean = NaN;
                end
                
                windbin = floor(WSmean/wind_interval)+1;
                if (windbin > num_wind_bins)
                    windbin = NaN;  %This will get excluded later
                end
                
                %Compute the radius of the sonde based on the current center
                %lat,lon here of the sonde are the LAUNCH coordinates
                [lat_center,lon_center,lat,lon,time_sonde] = get_track(F,[filedir files(i).name]);
                
                lat = lat-lat_center;
                lon = lon-lon_center;
                rearth = 6371; %[km]
                x = rearth*tand(lon);
                y = rearth*tand(lat);  %zoom in -- lon is in the "x" direction, lat is in the "y" direction
                rad = sqrt(x^2+y^2);
                
                %             %Compute the radial and tangential components:
                %             %Here is based on the lat/lon of LAUNCH:
                %             angle = -atan2(y,x);  %Make negative so we do the clockwise rotation
                %             rotmat = [cos(angle) -sin(angle); sin(angle) cos(angle)];
                %
                %             tmp_mat(1,:) = U;
                %             tmp_mat(2,:) = V;
                %
                %             radial_mat = rotmat*tmp_mat;
                %
                %             Ut = -radial_mat(1,:)'; %Tangential velocity -- CCW is positive direction
                %             Ur = radial_mat(2,:)'; %Radial direction -- outwards is positive
                
                
                %Compute the radial and tangential components:
                %Here is based on the lat/lon of CURRENT sonde reading
                for sonde_idx = 1:length(U)
                    lat_s = latvec(sonde_idx) - lat_center;
                    lon_s = lonvec(sonde_idx) - lon_center;
                    xs = rearth*tand(lon_s);
                    ys = rearth*tand(lat_s);
                    
                    
                    atan_tmp = atan2(ys,xs); %Location of lat/lon in polar angle
                    angle = atan_tmp + (atan_tmp<0)*2*pi;  %Need to correct for atan2 returning between -pi and pi
                    angle = angle - pi/2;  %Get the CCW rotation angle right according to the polar angle
                    rotmat = [cos(2*pi-angle) -sin(2*pi-angle); sin(2*pi-angle) cos(2*pi-angle)];
                    
                    tmp_mat = [U(sonde_idx); V(sonde_idx)];
                    radial_mat = rotmat*tmp_mat;
                    
                    Ut(sonde_idx) = -radial_mat(1);  %Negative to get positive Ut in the CCW direction
                    Ur(sonde_idx) = radial_mat(2);
                    
                end
                               
                clearvars tmp_mat radial_mat
                
                % Find RMW data closest to dropsonde in time, within 1 hr
                % for TCOBS & 6 hours for EBT
                if isTCOBS>0
                    k_rmw_time_hurr = find(rmw_time_hurr<=time_sonde+(1/48) & rmw_time_hurr>time_sonde-(1/48)); % Find TCOBS within an hour of the sonde
                elseif isEBT>0
                    k_rmw_time_hurr = find(rmw_time_hurr<=time_sonde+(1/8) & rmw_time_hurr>time_sonde-(1/8)); % Find EBT within 6 hours of the sonde
                end
                r_rmw_rad_hurr = -1;
                radbin = 0;
                rad_sonde_rms_ratio = -1;
                datevec(time_sonde);
                
                if (~isempty(k_rmw_time_hurr)) % representing whether time of hurricane has a corresponding RMW
                    r_rmw_rad_hurr  = rmw_rad_hurr(k_rmw_time_hurr); %radii of maximum wind speed of this dropsonde
                    v_rmw_speed_hurr = rmw_speed_hurr(k_rmw_time_hurr); %wind of maximum wind speed of this dropsonde
                    if (r_rmw_rad_hurr<100)  % from Dan C: if r_rmw_rad_hurr>100, very weak and disorganized storms
                        rad_sonde_rms_ratio = rad/r_rmw_rad_hurr; % dimentionless based on the radii ratio of RMW for different category hurricanes.
                        
                        radbin = floor(rad_sonde_rms_ratio/rad_interval)+1;
                        if (radbin > num_rad_bins)
                            radbin = 0; %Don't include if it falls outside of range
                        end
                        
                        %Guiquan's:
                        %                     if (rad_sonde_rms_ratio>0 && rad_sonde_rms_ratio<=1) %0-1RMW (eye)
                        %                        radbin=1;
                        %                     elseif (rad_sonde_rms_ratio>1 && rad_sonde_rms_ratio<=2)  %1-2RMW (away from wall)
                        %                        radbin=2;
                        %                     elseif (rad_sonde_rms_ratio>2 && rad_sonde_rms_ratio<=3)  %2-3RMW
                        %                        radbin=3;
                        %                     elseif (rad_sonde_rms_ratio>3)  %>3RMW
                        %                        radbin=4;
                        %                     end
                        
                    end % r_rmw_rad_hurr<100
                end % k_rmw_time_hurr from dropsonds has a corresponding RMW from Dan's data
                
                
                
                %Guiquan                theta_surf = theta(length(theta));
                %Guiquan                SST = get_sst([filedir files(i).name]);
                
                %Guiquan                thetatest = SST;
                %Guiquan                Dtheta = SST-theta_surf;
                
                %Guiquan                if (~isnan(SST))
                
                %Assign delta_theta bin:
                %Guiquan                thetabin = floor((thetatest-min_theta)/theta_interval)+1;
                %Guiquan                if (thetabin < 1); thetabin = 1; end
                %Guiquan                if (thetabin > num_theta_bins); thetabin = num_theta_bins; end
                
                %Guiquan                if (~isnan(radbin))
                %Guiquan                ssttosave{radbin,thetabin} = [ssttosave{radbin,thetabin} SST];
                %Guiquan                end
                
                
                %Now bin into the z-grid:
                if (radbin~=0 && ~isnan(radbin) && ~isnan(windbin))
                    
                    %Plot of all of the individual sondes in each radbin:
%                     figure(100)
%                     hold on
%                     subplot(floor(sqrt(num_rad_bins)),ceil(sqrt(num_rad_bins)),radbin)
%                     plot(WS,z,'-k')
%                     
%                     figure(101)
%                     hold on
%                     subplot(ceil(sqrt(num_wind_bins)),ceil(sqrt(num_wind_bins)),windbin)
%                     plot(WS,z,'-k')
                    
                    for j=1:length(z)
                        if (z(j)<max_height&&z(j)>min_height)
                            
%                             if (j > 1 && abs(WS(j)-WS(j-1)) > 10)
%                                 [filedir files(i).name]
%                                 figure(50)
%                                 hold on
%                                 plot(WS,z)
%                             end
                            
                            zidx = floor((z(j)-min_height)/height_interval)+1;
                            mean_U_profiles(zidx,radbin,windbin) = mean_U_profiles(zidx,radbin,windbin)+WS(j);
                            mean_Ur_profiles(zidx,radbin,windbin) = mean_Ur_profiles(zidx,radbin,windbin)+Ur(j);
                            mean_Ut_profiles(zidx,radbin,windbin) = mean_Ut_profiles(zidx,radbin,windbin)+Ut(j);
                            
                            mean_U_profiles_onestorm(zidx,radbin,windbin) = mean_U_profiles_onestorm(zidx,radbin,windbin)+WS(j);
                            mean_T_profiles(zidx,radbin,windbin) = mean_T_profiles(zidx,radbin,windbin)+T(j);
                            mean_RH_profiles(zidx,radbin,windbin) = mean_RH_profiles(zidx,radbin,windbin)+RH(j);
                            mean_zU_profiles(zidx,radbin,windbin) = mean_zU_profiles(zidx,radbin,windbin)+z(j);
                            mean_z_profiles(zidx,radbin,windbin) = mean_z_profiles(zidx,radbin,windbin)+z(j);
                            mean_p_profiles(zidx,radbin,windbin) = mean_p_profiles(zidx,radbin,windbin)+p(j);
                            mean_q_profiles(zidx,radbin,windbin) = mean_q_profiles(zidx,radbin,windbin)+q(j);
                            mean_k_profiles(zidx,radbin,windbin) = mean_k_profiles(zidx,radbin,windbin)+k(j);
                            mean_theta_profiles(zidx,radbin,windbin) = mean_theta_profiles(zidx,radbin,windbin)+theta(j);
                            
                            numvecU(zidx,radbin,windbin) = numvecU(zidx,radbin,windbin)+1;
                            numvecU_onestorm(zidx,radbin,windbin) = numvecU_onestorm(zidx,radbin,windbin)+1;
                            %Guiquan                        numvecT(zidx,radbin,thetabin) = numvecT(zidx,radbin,thetabin)+1;
                            
                            %Fill arrays with ALL data to take statistics:
                            all_U_profiles{zidx,radbin,windbin} = [all_U_profiles{zidx,radbin,windbin} WS(j)];
                            all_T_profiles{zidx,radbin,windbin} = [all_T_profiles{zidx,radbin,windbin} T(j)];
                            all_RH_profiles{zidx,radbin,windbin} = [all_RH_profiles{zidx,radbin,windbin} RH(j)];
                            all_z_profiles{zidx,radbin,windbin} = [all_z_profiles{zidx,radbin,windbin} z(j)];
                            all_zU_profiles{zidx,radbin,windbin} = [all_zU_profiles{zidx,radbin,windbin} z(j)];
                            all_p_profiles{zidx,radbin,windbin} = [all_p_profiles{zidx,radbin,windbin} p(j)];
                            all_q_profiles{zidx,radbin,windbin} = [all_q_profiles{zidx,radbin,windbin} q(j)];
                            all_k_profiles{zidx,radbin,windbin} = [all_k_profiles{zidx,radbin,windbin} k(j)];
                            all_theta_profiles{zidx,radbin,windbin} = [all_theta_profiles{zidx,radbin,windbin} theta(j)];
                        end
                    end %z-grid
                    
                    numsonde(radbin,windbin) = numsonde(radbin,windbin) + 1;
                    
                end
                
                
                sondexvec(count_sonde) = x/r_rmw_rad_hurr;
                sondeyvec(count_sonde) = y/r_rmw_rad_hurr;
                count_sonde = count_sonde + 1;
                
                %Guiquan                end %If ~(isnan(SST))
            end  %If length(z)>0 statement
        end %All the .frd files in one hurricane
       
        
    end % hurr in the rmw file?
end %Loop over hurricanes

fprintf('Total number of profiles used: %i\n',sum(sum(numsonde)));

mean_zU_profiles(:,:,:) = mean_zU_profiles(:,:,:)./numvecU(:,:,:);
mean_z_profiles(:,:,:)=mean_z_profiles(:,:,:)./numvecU(:,:,:);
mean_U_profiles(:,:,:) = mean_U_profiles(:,:,:)./numvecU(:,:,:);
mean_Ur_profiles(:,:,:) = mean_Ur_profiles(:,:,:)./numvecU(:,:,:);
mean_Ut_profiles(:,:,:) = mean_Ut_profiles(:,:,:)./numvecU(:,:,:);
mean_T_profiles(:,:,:)=mean_T_profiles(:,:,:)./numvecU(:,:,:);
mean_RH_profiles(:,:,:)=mean_RH_profiles(:,:,:)./numvecU(:,:,:);
mean_p_profiles(:,:,:)=mean_p_profiles(:,:,:)./numvecU(:,:,:);
mean_q_profiles(:,:,:)=mean_q_profiles(:,:,:)./numvecU(:,:,:);
mean_k_profiles(:,:,:)=mean_k_profiles(:,:,:)./numvecU(:,:,:);
mean_theta_profiles(:,:,:)=mean_theta_profiles(:,:,:)./numvecU(:,:,:);



%Plot the contours of mean velocity, as in Zhang et al. (2011):
[R,Z] = meshgrid(rplot,zplot);

% figure(2)
% hold on
% contourf(R,Z,mean_U_profiles,20,'edgecolor','none');
% caxis([0 60])
% axis([0 5 0 2000])
% %contour(R,Z,mean_U_profiles,0:5:60,'edgecolor','black')
% %[~, I] = max(mean_U_profiles,[],1);
% %plot(rplot,zplot(I),'-k','linewidth',3)
% title('WS')

figure(3)
hold on
plot(sondexvec,sondeyvec,'x')
xlabel('x')
ylabel('y')
axis([-max_rad max_rad -max_rad max_rad])

% figure(4)
% hold on
% contourf(R,Z,mean_Ur_profiles,20,'edgecolor','none');
% caxis([-20 15])
% axis([0 5 0 2000])
% title('Ur')
% 
% 
figure(5)
hold on
%contourf(R,Z,mean_Ut_profiles,20,'edgecolor','none');
caxis([0 60])
axis([0 5 0 2000])
title('Ut')

figure(6)
hold on
num_samples = sum(numsonde(:,:),2);
plot(rplot,num_samples,'-*k')

% figure(7)
% hold on
% contourf(R,Z,mean_theta_profiles,20,'edgecolor','none');
% caxis([302 320])
% axis([0 5 0 2000])
% title('T')
% 
% figure(8)
% hold on
% contourf(R,Z,mean_RH_profiles,20,'edgecolor','none');
% caxis([80 100])
% axis([0 5 0 2000])
% title('RH')

% figure(100)
% for ridx = 1:num_rad_bins
%     subplot(ceil(sqrt(num_rad_bins)),ceil(sqrt(num_rad_bins)),ridx)
%     title(['R/RMW = ' num2str(rplot(ridx))])
%     axis([-inf inf 0 2000])
% end
% sgtitle('Radius bins')
% 
% figure(101)
% for widx = 1:num_wind_bins
%     subplot(ceil(sqrt(num_wind_bins)),ceil(sqrt(num_wind_bins)),widx)
%     title(['WS = ' num2str(Uplot(widx))])
%     axis([-inf inf 0 2000])
% end
% sgtitle('Wind bins')

%% Saving Variables

%save('allStorms_all_sst_rad.mat','ssttosave');
save('allStorms_all_profile_data_rad.mat','all_U_profiles','all_theta_profiles','all_q_profiles','all_k_profiles','all_p_profiles');
save('allStorms_constants_rad.mat','rho','cpa','cpl','cpv','Lv','pbl_height','max_height',...
     'min_height','height_interval','num_z','max_rad','rad_interval',...
     'num_rad_bins','max_wind','wind_interval','num_wind_bins');


