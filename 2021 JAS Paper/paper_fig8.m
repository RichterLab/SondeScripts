close all
clear all
clc

%% Make azimuthally averaged slice of George's CM1 data
rad = squeeze(ncread('azim_avg.nc','xh'));
z = squeeze(ncread('azim_avg.nc','zh'));
u = squeeze(ncread('azim_avg.nc','uinterp'));
v = squeeze(ncread('azim_avg.nc','vinterp'));
wspd = sqrt(u.^2 + v.^2);

radct = 1;
for r=190:10:length(rad)
   umaxht(radct) = max(find(wspd(r,1:128) == nanmax(wspd(r,1:128))));
   radht(radct) = rad(r);
   radct = radct+1;
end
umaxhtmn = ceil(movmean(umaxht,9))

figure(7)
set(gcf,'position',[50 50 800 400])
pcolor(1e-3*rad,z,wspd')
shading flat
set(gca,'fontsize',14)
h = colorbar;
ylabel(h,'wind speed [m s^{-1}]','fontsize',14)
xlabel('radial distance [km]')
ylabel('height [m]')
caxis([0 90])
hold on
plot(1e-3*radht,z(umaxhtmn),'k--')