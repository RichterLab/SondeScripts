%%
clear all
close all
clc

%Load up the virtual sondes right up front:

% Filename containing trajectory data:
filename='./cm1out_pdata.nc';

time_traj=ncread(filename,'time');

tmp = ncread(filename,'z');
zfull(:,:) = tmp(:,:); %height
tmp = ncread(filename,'x');
xfull(:,:) = tmp(:,:); %make sure to check if the value is a 4D vector and then made sure same spot of 4Dmatrix are the defined ones (DONE!)
tmp = ncread(filename, 'y');
yfull(:,:) = tmp(:,:);
tmp = ncread(filename,'th');
thetafull(:,:) = tmp(:,:); % theta--take out T in other file
tmp = ncread(filename, 'prs');
presfull(:,:) = tmp(:,:); % pressure
tmp = ncread(filename,'u');
ufull(:,:)= tmp(:,:);
tmp = ncread(filename,'v');
vfull(:,:) = tmp(:,:);
wspdfull=sqrt(ufull.^2+vfull.^2); % windspeed!
tmp = ncread(filename, 'qv');
qfull(:,:) = tmp(:,:); %assume specific humidity (water vapor mixing ratio)


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

pbl_height = 10000;  %Height over which mean is taken (Powell uses 500m)
max_height = 10000;  %Height over which the fit is done
min_height = 10;   %Bottom height over which fit is done
height_interval = 10;
num_z = (max_height-min_height)/height_interval;

% Limits of wind speed bins
max_wind = 80;
wind_interval = 10;
num_wind_bins = max_wind/wind_interval;

%Limits of radius bins -- in normalized units (R/RMW)
max_rad = 5;
rad_interval = 0.5;
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



count_hurr=0;
count_sonde = 1;

count_hurr=count_hurr+1;
mean_U_profiles_onestorm = zeros(num_z,num_rad_bins,num_wind_bins);
numvecU_onestorm = zeros(num_z,num_rad_bins,num_wind_bins);



for i=1:1:numtraj % 2. loop every dropsonde
    
    
    %Assign the values - only take data where z,temp,p,WS,RH above 0
    ztmp = zfull(i,:); % [m]
    xtmp = xfull(i,:);
    ytmp = yfull(i,:);
    thetatmp = thetafull(i,:);
    ptmp = presfull(i,:); % [mb]
    WStmp = wspdfull(i,:); % [m/s]
    qtmp = qfull(i,:); % using for humidity
    utmp = ufull(i,:);
    vtmp = vfull(i,:);
    
    
    %Clean the data - if any of them have a negative signal, remove for
    %all profiles (only want a profile if all data is there)
    % This is the set I'll be plotting for now (all with full at end)
    z = ztmp(ztmp>0&thetatmp>0&ptmp>0&WStmp>0&qtmp>0);
    x = xtmp(ztmp>0&thetatmp>0&ptmp>0&WStmp>0&qtmp>0);
    y = ytmp(ztmp>0&thetatmp>0&ptmp>0&WStmp>0&qtmp>0);
    p = ptmp(ztmp>0&thetatmp>0&ptmp>0&WStmp>0&qtmp>0);
    q = qtmp(ztmp>0&thetatmp>0&ptmp>0&WStmp>0&qtmp>0);
    theta = thetatmp(ztmp>0&thetatmp>0&ptmp>0&WStmp>0&qtmp>0);
    WS = WStmp(ztmp>0&thetatmp>0&ptmp>0&WStmp>0&qtmp>0);
    u = utmp(ztmp>0&thetatmp>0&ptmp>0&WStmp>0&qtmp>0);
    v = vtmp(ztmp>0&thetatmp>0&ptmp>0&WStmp>0&qtmp>0);
    Ut = zeros(size(u));
    Ur = zeros(size(u));
    
    if (~isempty(z))  %Only do the computations if data is left after the clean
        
        
        %Now compute enthalpy k, as defined in Emanuel (1995):
        %k will have units of kJ/kg
        
        %Test using temperature:
        k = (cpa*(1-q) + cpl*q).*theta + Lv*q;
        
        %Compute mean velocity below pbl_height (Powell uses 500m):
        %WSmean = mean(WS(z<pbl_height));
        
        %Trying by using the lowest recorded wind speed rather than a PBL average
        if (z(end) < 100)  %Only do this if the lowest recording is beneath 100m
            WSmean = WS(end);
        else
            WSmean = NaN;
        end
        
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
        
        windbin = floor(double(WSmean)/wind_interval)+1;
        if (windbin > num_wind_bins)
            windbin = NaN;  %This will get excluded later
        end
        
        
        xcenter=0;
        ycenter=0;
        
        xdist=double(xfull(i,1))-xcenter;
        ydist=double(yfull(i,1))-ycenter;
        alpha=atan2(ydist,xdist);
        vt=cos(alpha).*vfull(i,1)-sin(alpha).*ufull(i,1);
        vr=cos(alpha).*ufull(i,1)+sin(alpha).*vfull(i,1);
        
        dist=sqrt(xdist.^2+ydist.^2); % r (radius from center)
        rad=dist/1000;
        
        r_rmw_rad_hurr  = 10; %rmw of the numerical simulation
        
        rad_sonde_rms_ratio = rad/r_rmw_rad_hurr; % dimentionless based on the radii ratio of RMW for different category hurricanes.
        
        radbin = floor(rad_sonde_rms_ratio/rad_interval)+1;
        if (radbin > num_rad_bins)
            radbin = 0; %Don't include if it falls outside of range
        end
        
        
        %Compute the radial and tangential components:
        %Here is based on the lat/lon of LAUNCH:
        %         angle = -atan2(ydist,xdist);  %Make negative so we do the clockwise rotation
        %         rotmat = [cos(angle) -sin(angle); sin(angle) cos(angle)];
        %
        %         tmp_mat(1,:) = u;
        %         tmp_mat(2,:) = v;
        %
        %         radial_mat = rotmat*tmp_mat;
        %
        %         Ut = -radial_mat(1,:)'; %Tangential velocity -- CCW is positive direction
        %         Ur = radial_mat(2,:)'; %Radial direction -- outwards is positive
        
        for sonde_idx = 1:length(u)
            xs = x(sonde_idx)/1000;
            ys = y(sonde_idx)/1000;
            
            
            atan_tmp = atan2(ys,xs); %Location of lat/lon in polar angle
            angle = atan_tmp + (atan_tmp<0)*2*pi;  %Need to correct for atan2 returning between -pi and pi
            angle = angle - pi/2;  %Get the CCW rotation angle right according to the polar angle
            rotmat = [cos(2*pi-angle) -sin(2*pi-angle); sin(2*pi-angle) cos(2*pi-angle)];
            
            tmp_mat = [u(sonde_idx); v(sonde_idx)];
            radial_mat = rotmat*tmp_mat;
            
            Ut(sonde_idx) = -radial_mat(1);  %Negative to get positive Ut in the CCW direction
            Ur(sonde_idx) = radial_mat(2);
            
        end
        
        
        clearvars tmp_mat
        
        
%         if (rad < 100)
%             figure(12)
%             clf
%             hold on
%             plot(u,z,'b')
%             plot(v,z,'g')
%             plot(Ur,z,'c')
%             plot(Ut,z,'r')
%             legend('U','V','Ur','Ut')
%             drawnow
%             
%             figure(13)
%             %clf
%             hold on
%             circlex = linspace(-11,11,100);
%             circletop = sqrt(11^2-circlex.^2);
%             circlebot = -sqrt(11^2-circlex.^2);
%             plot(circlex,circletop,'--k')
%             plot(circlex,circlebot,'--k')
%             qp=quiver(xdist/1000,ydist/1000,mean(u(z<200&z>50)),mean(v(z<200&z>50)));
%             %qp.MaxHeadSize = 0.8;
%             qp.Color = 'g';
%             axis([-100 100 -100 100])
%             drawnow
%             
%             pause
%         end
        
        
    end % isempty(z)
    
    
    
    %Now bin into the z-grid:
    if (radbin~=0 && ~isnan(radbin) && ~isnan(windbin))
        
        if (mod(i,10) == 0)
        %Plot of all of the individual sondes in each radbin:
        figure(100)
        hold on
        subplot(floor(sqrt(num_rad_bins)),ceil(sqrt(num_rad_bins)),radbin)
        plot(WS,z,'-k')
        
        figure(101)
        hold on
        subplot(ceil(sqrt(num_wind_bins)),ceil(sqrt(num_wind_bins)),windbin)
        plot(WS,z,'-k')
        end
        
        for j=1:length(z)
            if (z(j)<max_height&&z(j)>min_height)
                
                zidx = floor((z(j)-min_height)/height_interval)+1;
                mean_U_profiles(zidx,radbin,windbin) = mean_U_profiles(zidx,radbin,windbin)+WS(j);
                mean_Ur_profiles(zidx,radbin,windbin) = mean_Ur_profiles(zidx,radbin,windbin)+Ur(j);
                mean_Ut_profiles(zidx,radbin,windbin) = mean_Ut_profiles(zidx,radbin,windbin)+Ut(j);
                mean_U_profiles_onestorm(zidx,radbin,windbin) = mean_U_profiles_onestorm(zidx,radbin,windbin)+WS(j);
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
    
    
    sondexvec(count_sonde) = xdist/r_rmw_rad_hurr/1000;
    sondeyvec(count_sonde) = ydist/r_rmw_rad_hurr/1000;
    count_sonde = count_sonde + 1;
    
    %Guiquan                end %If ~(isnan(SST))
end  %all of the profiles


%Plot the contours of mean velocity, as in Zhang et al. (2011):
mean_U_profiles_onestorm(:,:,:) = mean_U_profiles_onestorm(:,:,:)./numvecU_onestorm(:,:,:);

% [R,Z] = meshgrid(rplot,zplot);
% figure(1)
% hold on
% contourf(R,Z,mean_U_profiles_onestorm,20,'edgecolor','none');
% caxis([0 80])
% axis([0 5 0 2000])
% colorbar
%contour(R,Z,mean_U_profiles_onestorm,0:5:60,'edgecolor','black')
%[maxspeed, I] = max(mean_U_profiles_onestorm,[],1);
%plot(rplot,zplot(I),'-k','linewidth',3)


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
% figure(5)
% hold on
% contourf(R,Z,mean_Ut_profiles,20,'edgecolor','none');
% caxis([0 60])
% axis([0 5 0 2000])
% title('Ut')

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

figure(100)
for ridx = 1:num_rad_bins
    subplot(floor(sqrt(num_rad_bins)),ceil(sqrt(num_rad_bins)),ridx)
    title(['R/RMW = ' num2str(rplot(ridx))])
    axis([-inf inf 0 2000])
end
sgtitle('Radius bins')

figure(101)
for widx = 1:num_wind_bins
    subplot(ceil(sqrt(num_wind_bins)),ceil(sqrt(num_wind_bins)),widx)
    title(['WS = ' num2str(Uplot(widx))])
    axis([-inf inf 0 2000])
end
sgtitle('Wind bins')

%% Saving Variables

%save('allStorms_all_sst_rad.mat','ssttosave');
save('allStorms_all_profile_data_rad.mat','all_U_profiles','all_theta_profiles','all_q_profiles','all_k_profiles','all_p_profiles');
save('allStorms_constants_rad.mat','rho','cpa','cpl','cpv','Lv','pbl_height','max_height',...
    'min_height','height_interval','num_z','max_rad','rad_interval',...
    'num_rad_bins','max_wind','wind_interval','num_wind_bins');


