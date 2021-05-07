close all
clear all
clc

load allProfiles_3km_vmax
all_U_profiles(find(all_U_profiles>110)) = NaN;
all_U_profiles(:,find(minZ>150)) = NaN;
all_U_profiles(:,find(stormVmax<50)) = NaN;

radRMW = radius_km./RMW;

badData = [1253 2175 4668 4828 5290 5296 5404 8221 6163 6287 6288 7373 9628 11434];
all_U_profiles(:,badData) = NaN; % Remove bad data, screened by Charlotte
all_z_profiles(:,badData) = NaN; % Remove bad data, screened by Charlotte
all_z_profiles(find(all_z_profiles==0))=NaN;

minRMW = [0:0.2:4.8];
maxRMW = [0.2:0.2:5];
zmin = [0:10:1990];
zmax = [10:10:2000];

for r = 1:length(minRMW)
    keep = find(radRMW > minRMW(r) & radRMW <= maxRMW(r));
    
    keepu = all_U_profiles(:,keep);
    keepu = reshape(keepu,[1 size(keepu,1)*size(keepu,2)]);
    keepz = all_z_profiles(:,keep);
    keepz = reshape(keepz,[1 size(keepz,1)*size(keepz,2)]);
    
    clear utmp
    utmp = nan(1,1);
    for z = 1:length(zmin)
       keepht = find(keepz > zmin(z) & keepz <= zmax(z));
       usave(r,z) = nanmean(keepu(keepht));
       zsave(r,z) = nanmean(keepz(keepht));
       
       clear keepht
    end
    clear keep keepu keepz
    r
end

% figure(1)
% pcolor(maxRMW,zmax,usave')
% shading flat,set(gca,'fontsize',12)
% shg
% h = colorbar;
% ylabel(h,'wind speed [m/s]','fontsize',12)
% ylabel('height [m]')
% xlabel('R/RMW')

for r=1:length(minRMW)
   umaxht(r) = find(usave(r,:)==nanmax(usave(r,:))) 
end

figure(2)
set(gcf,'position',[50 50 800 400])
contourf(maxRMW,zmax,usave',[0:1:80],'linecolor','none')
shading flat,set(gca,'fontsize',14)
shg
h = colorbar;
ylabel(h,'wind speed [m s^{-1}]','fontsize',12)
ylabel('height [m]')
xlabel('R/RMW')
hold on
umaxhtmn = movmean(umaxht,9)
plot(maxRMW,zmax(round(umaxhtmn)),'k--')
caxis([0 65])


 
 dataObs{2}.CData