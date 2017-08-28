
clear all
close all
clc

cpa = 1.005;  %kJ/kg-K
cpl = 4.19;   %kJ/kg-K
cpv = 1.87;   %kJ/kg-K
Lv = 2501;    %kJ/kg
rho = 1.2;   %kg/m^3
Rd = 0.28706; %kJ/kg-K

% Split data into bins for height, temp., windspeed, and relative humidity.
% Focus only on breaking up height bins first.

pbl_height = 500;  %Height over which mean is taken (Powell used 500m)
max_height = 100;  %Height over which the fit is done
min_height = 50;   %Bottom height over which fit is done
height_interval = 2; %Binning interval for height
num_z = (max_height-min_height)/height_interval;

%Limits of wind speed bins
max_wind = 80;
wind_interval = 80;
num_wind_bins = max_wind/wind_interval;

%Limits of the potential temperature bins:
min_theta = 301;
max_theta = 302;
theta_interval = 1;
num_theta_bins = (max_theta-min_theta)/theta_interval;

%Limits of the radius bins (km):
max_rad = 50e3;
min_rad = 0;
rad_interval = 5e3;
num_rad_bins = (max_rad-min_rad)/rad_interval;

%Making blank matrices to save spaces for data
numvecU = zeros(num_z,num_wind_bins,num_rad_bins);
numvecT = zeros(num_z,num_wind_bins,num_theta_bins,num_rad_bins);

all_U_profiles = cell(num_z,num_wind_bins,num_rad_bins);
all_T_profiles = cell(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
all_RH_profiles = cell(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
all_zU_profiles = cell(num_z,num_wind_bins,num_rad_bins);
all_z_profiles = cell(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
all_p_profiles = cell(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
all_q_profiles = cell(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
all_k_profiles = cell(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
all_theta_profiles = cell(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
ssttosave = cell(num_wind_bins,num_theta_bins,num_rad_bins);

% Filename containing trajectory data

%This is what we had used with Rachel:
%filename=['./cat5_050115.nc'];
%Nprofiles = 103041; % stopped shy of 103401 bc final nc file has errors


filename='./cm1out_pdata.nc';

time_traj=ncread(filename,'time');
zfull=ncread(filename,'z'); % height
xfull = ncread(filename,'x');
yfull = ncread(filename,'y');
thetafull=ncread(filename,'th'); % theta--take out T in other file
presfull=ncread(filename,'prs'); % pressure
ufull=ncread(filename,'u');
vfull=ncread(filename,'v');
wspdfull=sqrt(ufull.^2+vfull.^2); % windspeed! 
qfull=ncread(filename,'qv'); %assume specific humidity (water vapor mixing ratio)

numtraj=size(ufull,1);
numtimes=size(ufull,2);
for j=1:numtraj % Cleaning up data
    if min(zfull(j,:))<=1
        ind_end(j)=find(zfull(j,:)<=1e-4,1);
        zfull(j,ind_end(j):end)=-999;
        xfull(j,ind_end(j):end)=-999;
        yfull(j,ind_end(j):end)=-999;
        wspdfull(j,ind_end(j):end)=-999;
        thetafull(j,ind_end(j):end)=-999;
        presfull(j,ind_end(j):end)=-999;
        qfull(j,ind_end(j):end)=-999;
        time_to_sfc(j)=3+time_traj(end)-time_traj(1);
    end
end


for i=1:numtraj
    
    %Assign the values - only take data where z,temp,p,WS,RH above 0
    ztmp = zfull(i,:); % [m]
    xtmp = xfull(i,:);
    ytmp = yfull(i,:);
    thetatmp = thetafull(i,:);
    ptmp = presfull(i,:); % [mb]
    WStmp = wspdfull(i,:); % [m/s]
    qtmp = qfull(i,:); % using for humidity
    
    %Clean the data - if any of them have a negative signal, remove for
    %all profiles (only want a profile if all data is there)
    z = ztmp(ztmp>0&ptmp>0&WStmp>0&qtmp>0);
    x = xtmp(ztmp>0&ptmp>0&WStmp>0&qtmp>0);
    y = ytmp(ztmp>0&ptmp>0&WStmp>0&qtmp>0);
    theta = thetatmp(ztmp>0&ptmp>0&WStmp>0&qtmp>0);
    p = ptmp(ztmp>0&ptmp>0&WStmp>0&qtmp>0);
    WS = WStmp(ztmp>0&ptmp>0&WStmp>0&qtmp>0);
    q = qtmp(ztmp>0&ptmp>0&WStmp>0&qtmp>0);
    
    %Assuming the center of the storm is (x,y) = (0,0)
    rad = sqrt( x.^2 + y.^2);
    
    psurf = p(length(p));
    T=theta.*((p/psurf).^(Rd/cpa)); % can get T from theta here
    theta=theta.*((100000/psurf).^(-Rd/cpa));
    
    %A bit redundant but *pbl contains the profile from bottom to
    %pbl_height
    zpbl = z(z<pbl_height);   %% pbl=planetary boundary layer
    xpbl = x(z<pbl_height);
    ypbl = y(z<pbl_height);
    Tpbl = T(z<pbl_height);
    thetapbl = theta(z<pbl_height);
    ppbl = p(z<pbl_height);
    WSpbl = WS(z<pbl_height);
    qpbl = q(z<pbl_height);
    
    if (~isempty(z))  %Only do the computations if data is left after the clean

        %First compute enthalpy
        %k will have units of kJ/kg
        
        %Add a correction for adiabatic expansion with height
        %("potential enthalpy")
        %k_uncorrected = (cpa*(1-q) + cpl*q).*T + Lv*q;
        k = (cpa*(1-q) + cpl*q).*theta + Lv*q;

        %Compute mean velocity below pbl_height (Powell uses 500m):
        WSmean = mean(WSpbl);
        
        %Assign the velocity bin:
        WSbin = floor(WSmean/wind_interval)+1;
        if (WSbin > num_wind_bins); WSbin = num_wind_bins; end
        
        theta_surf = theta(length(theta));
        SST=301.15;
        thetatest = SST;
        
        %Assign delta_theta bin:
        thetabin = floor((thetatest-min_theta)/theta_interval)+1;
        if (thetabin < 1); thetabin = 1; end
        if (thetabin > num_theta_bins); thetabin = num_theta_bins; end
        
        %Assign the radius bin:
        %Use the initial radius to bin:
        radius = rad(1);
        radbin = floor(radius/rad_interval)+1;
        if (radbin > num_rad_bins); radbin = num_rad_bins; end
        
        if (~isnan(WSbin))
            ssttosave{WSbin,thetabin} = [ssttosave{WSbin,thetabin,radbin} SST];
        end
        
        %Now bin into the z-grid:
        for j=1:length(z)
            if (z(j)<max_height&&z(j)>min_height)
                if (~isnan(WSbin))
                    zidx = floor((z(j)-min_height)/height_interval)+1;
 
                    numvecU(zidx,WSbin,radbin) = numvecU(zidx,WSbin,radbin)+1;
                    numvecT(zidx,WSbin,thetabin,radbin) = numvecT(zidx,WSbin,thetabin,radbin)+1;
                    
                    %Fill arrays with ALL data to take statistics:
                    all_U_profiles{zidx,WSbin,radbin} = [all_U_profiles{zidx,WSbin,radbin} WS(j)];
                    all_T_profiles{zidx,WSbin,thetabin,radbin} = [all_T_profiles{zidx,WSbin,thetabin,radbin} T(j)];
                    all_zU_profiles{zidx,WSbin,radbin} = [all_z_profiles{zidx,WSbin,radbin} z(j)];
                    all_z_profiles{zidx,WSbin,thetabin,radbin} = [all_z_profiles{zidx,WSbin,thetabin,radbin} z(j)];
                    all_p_profiles{zidx,WSbin,thetabin,radbin} = [all_p_profiles{zidx,WSbin,thetabin,radbin} p(j)];
                    all_q_profiles{zidx,WSbin,thetabin,radbin} = [all_q_profiles{zidx,WSbin,thetabin,radbin} q(j)];
                    all_k_profiles{zidx,WSbin,thetabin,radbin} = [all_k_profiles{zidx,WSbin,thetabin,radbin} k(j)];
                    all_theta_profiles{zidx,WSbin,thetabin,radbin} = [all_theta_profiles{zidx,WSbin,thetabin,radbin} theta(j)];
                end
            end
        end %Number of vertical points (j-loop)
    end  %If length(z)>0 statement
end %All the profiles (i-loop)

% Saving Variables
save('all_sst.mat','ssttosave');
save('all_profile_data.mat','all_U_profiles','all_theta_profiles','all_q_profiles','all_k_profiles','all_p_profiles');
save('constants.mat','rho','cpa','cpl','cpv','Lv','Rd','pbl_height','max_height',...
    'min_height','height_interval','num_z','max_wind','wind_interval',...
    'num_wind_bins','min_theta','max_theta','theta_interval','num_theta_bins');
