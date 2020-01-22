%%
clear all
close all
clc

%%
%Split into 10m/s-wide bins from 0 to 80 m/s (each with its own mean profile)
%Split vertical into 5m-wide bins

cpa = 1.005;  %kJ/kg-K
cpl = 4.19;   %kJ/kg-K
cpv = 1.86;   %kJ/kg-K
Lv = 2260;    %kJ/kg
rho = 1.2;   %kg/m^3

% Split data into bins for height, temp., windspeed, and relative humidity.
% Focus only on breaking up height bins first.

pbl_height = 2000;  %Height over which mean is taken (Powell uses 500m)
max_height = 2000;  %Height over which the fit is done
min_height = 10;   %Bottom height over which fit is done
height_interval = 10;
num_z = (max_height-min_height)/height_interval;

% Limits of wind speed bins
max_wind = 80;                             
wind_interval = 80;
num_wind_bins = max_wind/wind_interval;

%Limits of radius bins -- in normalized units (R/RMW)
max_rad = 10;                             
rad_interval = 0.4;
num_rad_bins = max_rad/rad_interval;

%According to Guiquan, do it in more discrete intervals:
%num_rad_bins = 4; % in the eyewall or outside the eyewall
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

zplot = min_height+0.5*height_interval:height_interval:max_height-0.5*height_interval;
rplot = 0.5*rad_interval:rad_interval:max_rad-0.5*rad_interval;


%All of them that we have track data for:
% hurrvec = { ...
%     'Erika1997/' ...
%     'Bonnie1998/' 'Danielle1998/' 'Earl1998/' 'Georges1998/' 'Mitch1998/' ...
%     'Bret1999/' 'Dennis1999/' 'Dora1999/' 'Eugene1999/' 'Emily1999/' 'Floyd1999/' 'Irene1999/' ...
%     'Helene2006/' ...
%     'Felix2007/' 'Ingrid2007/' 'Karen2007/' ...
%     'Dolly2008/' 'Fay2008/' 'Gustav2008/' 'Hanna2008/' 'Ike2008/' 'Kyle2008/' 'Paloma2008/' ...
%     'Danny2009/' ...
%     'Earl2010/' 'Karl2010/' 'Alex2010/' ...
%     'Dora2011/' 'Irene2011/' 'Rina2011/' ...
%     'Leslie2012/' 'Sandy2012/' ...
%     'Gabrielle2013/' 'Ingrid2013/' 'Karen2013/' ...
%     'Arthur2014/' 'Bertha2014/' 'Cristobal2014/' 'Edouard2014/' 'Gonzalo2014/' ...
%     'Danny2015/' 'Erika2015/' 'Guillermo2015/' 'Hilda2015/' 'Joaquin2015/' 'Oho2015/' ...
%     'Hermine2016/' 'Karl2016/' 'Matthew2016/' 'TropDep2016/' ...
%     'Franklin2017/' 'Harvey2017/' 'Irma2017/' 'Jose2017/' 'Maria2017/' 'Nate2017/'};
 
%Without basically empty ones (as determined by composite wind height-radius plots):
hurrvec = { ...
    'Erika1997/' ...
    'Bonnie1998/' 'Danielle1998/' 'Earl1998/' 'Georges1998/' 'Mitch1998/' ...
    'Bret1999/' 'Dennis1999/' 'Floyd1999/' ...
    'Helene2006/' ...
    'Felix2007/' 'Ingrid2007/' ...
    'Dolly2008/' 'Fay2008/' 'Gustav2008/' 'Ike2008/' 'Paloma2008/' ...
    'Earl2010/' 'Karl2010/' ...
    'Irene2011/' 'Rina2011/' ...
    'Leslie2012/' 'Sandy2012/' ...
    'Ingrid2013/' 'Karen2013/' ...
    'Arthur2014/' 'Bertha2014/' 'Cristobal2014/' 'Edouard2014/' 'Gonzalo2014/' ...
    'Danny2015/' 'Erika2015/' 'Joaquin2015/' ...
    'Hermine2016/' 'Karl2016/' 'Matthew2016/' ...
    'Franklin2017/' 'Harvey2017/' 'Irma2017/' 'Jose2017/' 'Maria2017/' 'Nate2017/'};

%Testing:
%hurrvec = {'Franklin2017/' 'Harvey2017/' 'Irma2017/' 'Jose2017/' 'Maria2017/' 'Nate2017/' };

%%% 1. RMW from Dan Chavas, every 6 hours%%%%%
load('./RWS/ebtrk_atlc_1988_2017.mat','Year_EBT','Month_EBT','Day_EBT','HourUTC_EBT','StormName_EBT','rmkm_EBT','Vmms_EBT');

%%% 2. RMW from Jonathan Vigh, every 1 hour%%%%%
%%% Attention: The datatype is different with 'ebtrk_atlc_1988_2017.mat'%%%
% load('./RWS/TCOBS_0.40_20160119.nc.mat');

cach = num2str(Year_EBT); %extract year added to rmw_name
rmw_name=strcat(StormName_EBT,cach,'/');
rmw_time=Month_EBT*10000+Day_EBT*24+HourUTC_EBT; %Form a 'number' to represent the rmw_time 
clearvars cach

profilecount(1:num_wind_bins,1:num_rad_bins,length(hurrvec))=0;

count_hurr=0;
count_sonde = 1;
for hurr = hurrvec % 1. loop every hurricane
    count_hurr=count_hurr+1;
    mean_U_profiles_onestorm = zeros(num_z,num_rad_bins,num_wind_bins);
    numvecU_onestorm = zeros(num_z,num_rad_bins,num_wind_bins);
    
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

        %A bit redundant but *pbl contains the profile from bottom to pbl_height (500m,used by Powell 2003)
        zpbl = z(z<pbl_height);   
        Tpbl = T(z<pbl_height);
        ppbl = p(z<pbl_height);
        WSpbl = WS(z<pbl_height);
        RHpbl = RH(z<pbl_height);

        if (~isempty(z))  %Only do the computations if data is left after the clean

            %Compute potential temperature:
            psurf = p(length(p));
            theta = T.*(psurf./p).^(2/7);
            thetapbl = Tpbl.*(psurf./ppbl).^(2/7);

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
            WSmean = mean(WSpbl);
            
            windbin = floor(WSmean/wind_interval)+1;
            if (windbin > num_wind_bins)
                windbin = num_wind_bins;
            end

            %Compute the radius of the sonde based on the current center
            [lat_center,lon_center,lat,lon,time_sonde] = get_track(hurr{1},[filedir files(i).name]);
            lat = lat-lat_center; lon = lon-lon_center;
            rearth = 6371; %[km]
            x = rearth*tand(lat); y = rearth*tand(lon);
            rad = sqrt(x^2+y^2);
    
            % position of dropsonde, because RMW is every 6 bours, so if the
            % k_rmw_time_hurr (time_sonde-3, time_sonde+3), then
            % r_rmw_rad_hurr = RMW(time sonde), without interpolation
            k_rmw_time_hurr = find(rmw_time_hurr<=time_sonde+3 & rmw_time_hurr>time_sonde-3); 
            r_rmw_rad_hurr = -1; 
            radbin = 0; 
            rad_sonde_rms_ratio = -1;
            
            if (~isempty(k_rmw_time_hurr)) % representing whether time of hurricane has a corresponding RMW 
                r_rmw_rad_hurr  = rmw_rad_hurr(k_rmw_time_hurr); %radii of maximum wind speed of this dropsonde
                v_rmw_speed_hurr = rmw_speed_hurr(k_rmw_time_hurr); %wind of maximum wind speed of this dropsonde
                if (r_rmw_rad_hurr<100)  % from Dan C: if r_rmw_rad_hurr>100, very weak and disorganized storms
                    rad_sonde_rms_ratio = rad/r_rmw_rad_hurr; % dimentionless based on the radii ratio of RMW for different category hurricanes.
                  
                    radbin = floor(rad_sonde_rms_ratio/rad_interval)+1;
                    if (radbin > num_rad_bins)
                        radbin = num_rad_bins;
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
            for j=1:length(z)
                if (z(j)<max_height&&z(j)>min_height)

                    zidx = floor((z(j)-min_height)/height_interval)+1;
                    mean_U_profiles(zidx,radbin,windbin) = mean_U_profiles(zidx,radbin,windbin)+WS(j);
                    mean_U_profiles_onestorm(zidx,radbin,windbin) = mean_U_profiles_onestorm(zidx,radbin,windbin)+WS(j);
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
                    numvecU_onestorm(zidx,radbin,windbin) = numvecU_onestorm(zidx,radbin,windbin)+1;
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
            end %z-grid
            end 

            if (radbin~=0 && ~isnan(radbin) && ~isnan(windbin)) 
                profilecount(windbin,radbin,count_hurr) = profilecount(windbin,radbin,count_hurr) + 1;
            end
            
            sondexvec(count_sonde) = x/r_rmw_rad_hurr;
            sondeyvec(count_sonde) = y/r_rmw_rad_hurr;
            count_sonde = count_sonde + 1;

%Guiquan                end %If ~(isnan(SST))
        end  %If length(z)>0 statement 
    end %All the .frd files in one hurricane
    
    %Plot the contours of mean velocity, as in Zhang et al. (2011):
    mean_U_profiles_onestorm(:,:,:) = mean_U_profiles_onestorm(:,:,:)./numvecU_onestorm(:,:,:);

    [R,Z] = meshgrid(rplot,zplot);
    figure(1)
    subplot(floor(sqrt(length(hurrvec))),ceil(sqrt(length(hurrvec))),count_hurr)
    hold on
    contourf(R,Z,mean_U_profiles_onestorm,0:2:60,'edgecolor','none');
    axis([0 5 0 2000])
    colorbar
    %contour(R,Z,mean_U_profiles_onestorm,0:5:60,'edgecolor','black')
    %[maxspeed, I] = max(mean_U_profiles_onestorm,[],1);
    %plot(rplot,zplot(I),'-k','linewidth',3)
    title(hurr)
    
 end % hurr in the rmw file?
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



%Plot the contours of mean velocity, as in Zhang et al. (2011):
[R,Z] = meshgrid(rplot,zplot);

figure(2)
hold on
contourf(R,Z,mean_U_profiles,0:0.1:60,'edgecolor','none');
axis([0 5 0 2000])
contour(R,Z,mean_U_profiles,0:5:60,'edgecolor','black')
[maxspeed, I] = max(mean_U_profiles,[],1);
plot(rplot,zplot(I),'-k','linewidth',3)

figure(3)
hold on
plot(sondexvec,sondeyvec,'x')
xlabel('x')
ylabel('y')
axis([-5 5 -5 5])

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

