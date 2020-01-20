%%
clear all
close all
clc

%%
%Split into 10m/s-wide bins from 0 to 80 m/s (each with its own mean profile)
%Split vertical into 5m-wide bins

%Guiquan cpa = 1.005;  %kJ/kg-K
%Guiquan cpl = 4.19;   %kJ/kg-K
%Guiquan cpv = 1.86;   %kJ/kg-K
%Guiquan Lv = 2260;    %kJ/kg
rho = 1.2;   %kg/m^3

% Split data into bins for height, temp., windspeed, and relative humidity.
% Focus only on breaking up height bins first.

pbl_height = 100;  %Height over which mean is taken (Powell uses 500m)
max_height = 2000;  %Height over which the fit is done
min_height = 10;   %Bottom height over which fit is done
height_interval = 5;
num_z = (max_height-min_height)/height_interval;

% Limits of wind speed bins
max_wind = 80;                             
wind_interval = 10;
num_wind_bins = max_wind/wind_interval;

%Limits of radius bins
max_rad = 300;                             
rad_interval = 60;
%num_rad_bins = max_rad/rad_interval;
num_rad_bins = 4; % in the eyewall or outside the eyewall
%0-0.5; 0.1-1; 1-2; 2-3; 3-4; >4
%0-1; 1-2; 2-3; >3

%Limits of the potential temperature bins:
%Guiquan min_theta = 295;                           
%Guiquan max_theta = 305;
%Guiquan theta_interval = 1;
%Guiquan num_theta_bins = (max_theta-min_theta)/theta_interval;

% Making blank matrices to save spaces for data
% [height_bin, radii_bin, wind_bin]
mean_U_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
numvecU = zeros(num_z,num_rad_bins,num_wind_bins);
%Guiquan mean_T_profiles = zeros(num_z,num_rad_bins,num_theta_bins);
%Guiquan numvecT = zeros(num_z,num_rad_bins,num_theta_bins);
%Guiquan mean_RH_profiles = zeros(num_z,num_rad_bins,num_theta_bins);
mean_zU_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
%Guiquan mean_z_profiles = zeros(num_z,num_rad_bins,num_theta_bins);
%Guiquan mean_p_profiles = zeros(num_z,num_rad_bins,num_theta_bins);
%Guiquan mean_q_profiles = zeros(num_z,num_rad_bins,num_theta_bins);
%Guiquan mean_k_profiles = zeros(num_z,num_rad_bins,num_theta_bins);
%Guiquan mean_theta_profiles = zeros(num_z,num_rad_bins,num_theta_bins);

all_U_profiles = cell(num_z,num_rad_bins,num_wind_bins);
%Guiquan all_T_profiles = cell(num_z,num_rad_bins,num_theta_bins);
%Guiquan all_RH_profiles = cell(num_z,num_rad_bins,num_theta_bins);
all_zU_profiles = cell(num_z,num_rad_bins,num_wind_bins);
%Guiquan all_z_profiles = cell(num_z,num_rad_bins,num_theta_bins);
%Guiquan all_p_profiles = cell(num_z,num_rad_bins,num_theta_bins);
%Guiquan all_q_profiles = cell(num_z,num_rad_bins,num_theta_bins);
%Guiquan all_k_profiles = cell(num_z,num_rad_bins,num_theta_bins);
%Guiquan all_theta_profiles = cell(num_z,num_rad_bins,num_theta_bins);
%Guiquan ssttosave = cell(num_rad_bins,num_theta_bins);

% hurrvec = {'Earl2010/'}; %'Bonnie1998/' 'Danielle1998/' 'Danny1997/'
% 'Dora1999/' 'Earl1998/' 'Erika1997/' 'Georges1998/' 'Guillermo1997/' 'Mitch1998/' 
% 'TD22010/' 'TD51997/' All these other ones b/c they don't have track
% files

% hurrvec = {'Earl2010/' 'Karl2010/' 'Alex2010/' ...
%     'Bret1999/' 'Danny2009/' 'Dennis1999/' 'Dolly2008/' 'Dora2011/'  ...
%     'Emily1999/' 'Eugene1999/' 'Fay2008/' 'Felix2007/' 'Floyd1999/'  ...
%     'Gustav2008/' 'Hanna2008/' 'Helene2006/' 'Ike2008/' 'Ingrid2007/' 'Irene1999/' 'Irene2011/' ...
%     'Karen2007/' 'Kyle2008/' 'Leslie2012/' 'Paloma2008/' 'Rina2011/' 'Sandy2012/' ...
%     'Gabrielle2013/' 'Ingrid2013/' 'Karen2013/' ...
%     'Arthur2014/' 'Bertha2014/' 'Cristobal2014/' 'Edouard2014/' 'Gonzalo2014/'...
%     'Danny2015/' 'Erika2015/' 'Guillermo2015/' 'Hilda2015/' 'Joaquin2015/' 'Oho2015/' ...
%     'Hermine2016/' 'Karl2016/' 'Matthew2016/' 'TropDep2016/' ...
%     'Franklin2017/' 'Harvey2017/' 'Irma2017/' 'Jose2017/' 'Maria2017/' 'Nate2017/'}; %updated dropsondes

hurrvec = {'Earl2010/' 'Maria2017/'};

%%% 1. RMW from Dan Chavas, every 6 hours%%%%%
load('./RWS/ebtrk_atlc_1988_2017.mat','Year_EBT','Month_EBT','Day_EBT','HourUTC_EBT','StormName_EBT','rmkm_EBT','Vmms_EBT');
%%% 2. RMW from Dan Chavas, every 1 hour%%%%%
%%% Attention: The datatype is different with 'ebtrk_atlc_1988_2017.mat'%%%
% load('./RWS/TCOBS_0.40_20160119.nc.mat');
cach = num2str(Year_EBT); %extract year added to rmw_name
rmw_name=strcat(StormName_EBT,cach,'/');
rmw_time=Month_EBT*10000+Day_EBT*24+HourUTC_EBT; %Form a 'number' to represent for the rmw_time 
clearvars cach
%%
profilecount(1:num_wind_bins,1:num_rad_bins,length(hurrvec))=0;
% profilecount(1:8,1:4,1:48) = 0;
count_hurr=0;
for hurr = hurrvec % 1. loop every hurricane
    count_hurr=count_hurr+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cach = strmatch(upper(hurr), rmw_name); % RMW location of specific hurricane
 if (~isnan(cach)) % ensure hurr name of dropsondes corresponds to a rmw_name from the track data
    clearvars rmw_time_hurr rmw_rad_hurr 
    cach1=rmw_time(cach(1):cach(end)); % RMW time of specific hurricane
    cach2=rmkm_EBT(cach(1):cach(end)); % RMW radius of specific hurricane
    cach3=Vmms_EBT(cach(1):cach(end)); % RMW wind speed of specific hurricane
    rmw_time_hurr = cach1(cach1>0&cach2>0);%remove nan
    rmw_rad_hurr  = cach2(cach1>0&cach2>0);%remove nan
    rmw_speed_hurr= cach3(cach1>0&cach2>0);%remove nan
    clearvars cach1 cach2 cach3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    filedirtmp = strcat('./all_storms/',hurr);
    filedir = filedirtmp{1};
    files = dir([filedir '*.frd']);
    for i=1:length(files) % 2. loop every dropsonde

        fprintf('Opening file: %s\n',[filedir files(i).name]);
        datstruct = importdata([filedir files(i).name],' ',21);
        dat = datstruct.data;

        %Assign the values - only take data where z,temp,p,WS,RH above 0
        ztmp = dat(:,6); % [m]
        Ttmp = dat(:,4); % [C]
        ptmp = dat(:,3); % [mb]
        WStmp = dat(:,8); % [m/s]
        RHtmp = dat(:,5); % [%]

        %Clean the data - if any of them have a negative signal, remove for
        %all profiles (only want a profile if all data is there)
        % This is the set I'll be plotting for now (all with full at end)
        z = ztmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0);
        T = Ttmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0) + 273;
        p = ptmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0);
        WS = WStmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0);
        RH = RHtmp(ztmp>0&Ttmp>0&ptmp>0&WStmp>0&RHtmp>0);

        %A bit redundant but *pbl contains the profile from bottom to
        %pbl_height (500m,used by Powell 2003)
        zpbl = z(z<pbl_height & z>10);   %% pbl=planetary boundary layer
        Tpbl = T(z<pbl_height & z>10);
        ppbl = p(z<pbl_height & z>10);
        WSpbl = WS(z<pbl_height & z>10);
        RHpbl = RH(z<pbl_height & z>10);

        if (~isempty(z))  %Only do the computations if data is left after the clean

            %Compute potential temperature:
%Guiquan                psurf = p(length(p));
%Guiquan                theta = T.*(psurf./p).^(2/7);
%Guiquan                thetapbl = Tpbl.*(psurf./ppbl).^(2/7);

            %Compute the specific humidity q:

            %Need saturation vapor pressure:
            %Formula from Dingman's hydrology book, eq. D-7 ("Magnus relation")
            %Needs T in celsius, gives estar in mb (hPa)
%Guiquan                estar = 6.11*exp(17.3*(T-273)./((T-273)+237.3));
%Guiquan                e = (RH/100).*estar;

            %p is in mb as well, so can divide:
%Guiquan                eoverp = e./p;
%Guiquan                q = 0.622*eoverp./(1-0.378*eoverp);

            %Now compute enthalpy k, as defined in Emanuel (1995):
            %k will have units of kJ/kg

            %Test using temperature:
%Guiquan                k_uncorrected = (cpa*(1-q) + cpl*q).*theta + Lv*q;

            %Add a correction for adiabatic expansion with height
            %("potential enthalpy")
%Guiquan                k = k_uncorrected;   %Realized I hadn't done this correctly at first - just need to use "theta" not "T" in k_uncorrected

            %Compute mean velocity below pbl_height (Powell uses 500m):
            WSmean = mean(WSpbl);

%               %Assign the velocity bin: 1=0-5m/s, 2=5-10m/s,etc.
%               WSbin = floor(WSmean/wind_interval)+1;
%               if (WSbin > num_wind_bins); WSbin = num_wind_bins; end

            % Assign the radius bin: 1=1-10m, 2=10-20m, etc.

            [lat_center,lon_center,lat,lon,time_sonde] = get_track(hurr{1},[filedir files(i).name]);
            lat = lat-lat_center; lon = lon-lon_center;
            rearth = 6371; %[km]
            x = rearth*tand(lat); y = rearth*tand(lon);
            rad = sqrt(x^2+y^2);
            
            k_rmw_time_hurr=-1;
            % position of dropsonde, because RMW is every 6 bours, so if the
            % k_rmw_time_hurr (time_sonde-3, time_sonde+3), then
            % r_rmw_rad_hurr = RMW(time sonde), without interpolation
            k_rmw_time_hurr=find(rmw_time_hurr<=time_sonde+3 & rmw_time_hurr>time_sonde-3); 
            r_rmw_rad_hurr=-1; 
            radbin =0; 
            rad_sonde_rms_ratio=-1;
            if (k_rmw_time_hurr~=-1) % representing for this time of hurricane has a corresponding RMW 
                r_rmw_rad_hurr  =rmw_rad_hurr(k_rmw_time_hurr); %radii of maximum wind speed of this dropsonde
                v_rmw_speed_hurr=rmw_speed_hurr(k_rmw_time_hurr); %wind of maximum wind speed of this dropsonde
                if (r_rmw_rad_hurr<100)% & v_rmw_speed_hurr>30) % from Dan, if r_rmw_rad_hurr>100, very weak and disorganized storms
                    rad_sonde_rms_ratio=rad/r_rmw_rad_hurr; % dimentionless based on the radii ratio of RMW for different category hurricanes.
                    if (rad_sonde_rms_ratio>0 && rad_sonde_rms_ratio<=1) %0-1RMW (eye)
                       radbin=1;
                    elseif (rad_sonde_rms_ratio>1 && rad_sonde_rms_ratio<=2)  %1-2RMW (away from wall)
                       radbin=2;
                    elseif (rad_sonde_rms_ratio>2 && rad_sonde_rms_ratio<=3)  %2-3RMW
                       radbin=3;
                    elseif (rad_sonde_rms_ratio>3)  %>3RMW
                       radbin=4;
                    end
                end % r_rmw_rad_hurr<100
            end % k_rmw_time_hurr from dropsonds has a corresponding RMW from Dan's data 
            
                      
            windbin = floor(WSmean/wind_interval)+1;
            if (windbin > num_wind_bins); windbin = num_wind_bins; end

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
            for j=1:length(z)
                if (z(j)<max_height&&z(j)>min_height)
                %if (~isnan(radbin) && ~isnan(windbin))
                if (radbin~=0 && ~isnan(windbin))
                    zidx = floor((z(j)-min_height)/height_interval)+1;
                    mean_U_profiles(zidx,radbin,windbin) = mean_U_profiles(zidx,radbin,windbin)+WS(j);
%Guiquan                        mean_T_profiles(zidx,radbin,thetabin) = mean_T_profiles(zidx,radbin,thetabin)+T(j);
%Guiquan                        mean_RH_profiles(zidx,radbin,thetabin) = mean_RH_profiles(zidx,radbin,thetabin)+RH(j);
                    mean_zU_profiles(zidx,radbin,windbin) = mean_zU_profiles(zidx,radbin,windbin)+z(j);
%Guiquan                        mean_z_profiles(zidx,radbin,thetabin) = mean_z_profiles(zidx,radbin,thetabin)+z(j);
%Guiquan                        mean_p_profiles(zidx,radbin,thetabin) = mean_p_profiles(zidx,radbin,thetabin)+p(j);
%Guiquan                        mean_q_profiles(zidx,radbin,thetabin) = mean_q_profiles(zidx,radbin,thetabin)+q(j);
%Guiquan                        mean_k_profiles(zidx,radbin,thetabin) = mean_k_profiles(zidx,radbin,thetabin)+k(j);
%Guiquan                        mean_theta_profiles(zidx,radbin,thetabin) = mean_theta_profiles(zidx,radbin,thetabin)+theta(j);
                    %radbin
                    numvecU(zidx,radbin,windbin) = numvecU(zidx,radbin,windbin)+1;
%Guiquan                        numvecT(zidx,radbin,thetabin) = numvecT(zidx,radbin,thetabin)+1;

                    %Fill arrays with ALL data to take statistics:
                    all_U_profiles{zidx,radbin,windbin} = [all_U_profiles{zidx,radbin,windbin} WS(j)];
%Guiquan                        all_T_profiles{zidx,radbin,thetabin} = [all_T_profiles{zidx,radbin,thetabin} T(j)];
%Guiquan                        all_RH_profiles{zidx,radbin,thetabin} = [all_RH_profiles{zidx,radbin,thetabin} RH(j)];
%Guiquan                        all_z_profiles{zidx,radbin,thetabin} = [all_z_profiles{zidx,radbin,thetabin} z(j)];
                    all_zU_profiles{zidx,radbin,windbin} = [all_zU_profiles{zidx,radbin,windbin} z(j)];
%Guiquan                        all_p_profiles{zidx,radbin,thetabin} = [all_p_profiles{zidx,radbin,thetabin} p(j)];
%Guiquan                        all_q_profiles{zidx,radbin,thetabin} = [all_q_profiles{zidx,radbin,thetabin} q(j)];
%Guiquan                        all_k_profiles{zidx,radbin,thetabin} = [all_k_profiles{zidx,radbin,thetabin} k(j)];
%Guiquan                        all_theta_profiles{zidx,radbin,thetabin} = [all_theta_profiles{zidx,radbin,thetabin} theta(j)];
                end
                end
            end %z-grid

            if (~isnan(windbin) && radbin >0) 
                profilecount(windbin,radbin,count_hurr) = profilecount(windbin,radbin,count_hurr) + 1;
            end

%Guiquan                end %If ~(isnan(SST))
        end  %If length(z)>0 statement 
    end %All the .frd files in one hurricane
 end % hurr in the rmw file
end %Loop over hurricanes

fprintf('Total number of profiles used: %i\n',sum(sum(sum(profilecount))));

mean_zU_profiles(:,:,:)=mean_zU_profiles(:,:,:)./numvecU(:,:,:);
%Guiquan    mean_z_profiles(:,:,:)=mean_z_profiles(:,:,:)./numvecT(:,:,:);
mean_U_profiles(:,:,:)=mean_U_profiles(:,:,:)./numvecU(:,:,:);
%Guiquan    mean_T_profiles(:,:,:)=mean_T_profiles(:,:,:)./numvecT(:,:,:);
%Guiquan    mean_RH_profiles(:,:,:)=mean_RH_profiles(:,:,:)./numvecT(:,:,:);
%Guiquan    mean_p_profiles(:,:,:)=mean_p_profiles(:,:,:)./numvecT(:,:,:);
%Guiquan    mean_q_profiles(:,:,:)=mean_q_profiles(:,:,:)./numvecT(:,:,:);
%Guiquan    mean_k_profiles(:,:,:)=mean_k_profiles(:,:,:)./numvecT(:,:,:);
%Guiquan    mean_theta_profiles(:,:,:)=mean_theta_profiles(:,:,:)./numvecT(:,:,:);
    
%% Saving Variables

%Guiquan save('allStorms_all_sst_rad.mat','ssttosave');
%Guiquan save('allStorms_all_profile_data_rad.mat','all_U_profiles','all_theta_profiles','all_q_profiles','all_k_profiles','all_p_profiles');
%Guiquan save('allStorms_constants_rad.mat','rho' ,'cpa','cpl','cpv','Lv','pbl_height','max_height',...
%Guiquan     'min_height','height_interval','num_z','max_rad','rad_interval',...
%Guiquan     'num_rad_bins','min_theta','max_theta','theta_interval','num_theta_bins');

save('allStorms_all_profile_data_rad_MBL100.mat','all_U_profiles');
save('allStorms_constants_rad_MBL100.mat','rho' ,'pbl_height','max_height',...
    'min_height','height_interval','num_z','max_rad','rad_interval',...
    'num_rad_bins','profilecount');
%% count the hurricane name, year, how many profiles from dropsonds
% in these frofiles, how many are in wind speed 0-10, 10-20, 20-30 ...................
fileID = fopen('name_year_MBL_h100.txt','w');
    fprintf(fileID,'%10s %6s %10s %6s %6s %6s %6s %6s %6s %6s %6s \n','name','year','Profiles','0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80');
for i=1:1:length(hurrvec)
    hurr=hurrvec{i};
    hurr_length=length(hurr);
    hurr_year=hurr(end-4:end-1);
    hurr_name=hurr(1:end-5);
    fprintf(fileID,'%10s %6s %10d %6d %6d %6d %6d %6d %6d %6d %6d \n',hurr_name,hurr_year,sum(sum(profilecount(:,:,i),2),1), sum(profilecount(:,:,i),2));
end
    fprintf(fileID,'%10s %6s %10d %6d %6d %6d %6d %6d %6d %6d %6d \n','total','--',sum(sum(sum(profilecount(:,:,:),2),1),3), sum(sum(profilecount(:,:,:),2),3));
fclose(fileID);
%% count how many profiles of [wind-speed, R/RMW]
fileID = fopen('RMW_MBL_h100.txt','w');
fprintf(fileID,'%10s %6s %6s %6s %6s %6s %6s %6s %6s \n','R/RMW','0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80');
fprintf(fileID,'%10s %6d %6d %6d %6d %6d %6d %6d %6d \n','0-1', sum(profilecount(:,1,:),3));
fprintf(fileID,'%10s %6d %6d %6d %6d %6d %6d %6d %6d \n','1-2', sum(profilecount(:,2,:),3));
fprintf(fileID,'%10s %6d %6d %6d %6d %6d %6d %6d %6d \n','2-3', sum(profilecount(:,3,:),3));
fprintf(fileID,'%10s %6d %6d %6d %6d %6d %6d %6d %6d \n','>3' , sum(profilecount(:,4,:),3));
fprintf(fileID,'%10s %6d %6d %6d %6d %6d %6d %6d %6d \n','Profiles', sum(sum(profilecount(:,1:4,:),3),2));
fclose(fileID);
%%
% zplot = min_height+0.5*height_interval:height_interval:max_height-0.5*height_interval;
% figure
% for i=1:1:16
%     a=squeeze(mean_U_profiles(:,2,i));
%     plot(a,zplot,'o-');
%     hold on
% end
