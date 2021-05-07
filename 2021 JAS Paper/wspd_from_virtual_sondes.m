%% This is needed to produce the dashed line on Fig. 8

close all
clear all
clc

load ('/Volumes/Elements/Work Projects/Hurricanes/virtualSondes/virtual_sonde_data_31m.mat')

RMW = 12e3;

for t=1:length(zdata)
    ztmp = zdata{t,:};
    %    wstmp = data.WSdata{t,:};
    xtmp = xdata{t,:};
    ytmp = ydata{t,:};
    zkeep = find(ztmp<=2000);
    %    meanWS500(t) = nanmean(wstmp(zkeep));
    %    maxWS(t) = nanmax(wstmp(zkeep));
    %    radius_top(t) = sqrt(xtmp(1).^2+ytmp(1).^2);
    rad(t) = nanmean(sqrt(xtmp.^2+ytmp.^2));
    clear wstmp ztmp xtmp ytmp
    t
end

radRMW = rad./RMW;

minRMW = [0:0.2:4.8];
maxRMW = [0.2:0.2:5];
zmin = [0:10:1990];
zmax = [10:10:2000];

wspd = zeros(25,200);
nwspd = zeros(25,200);
for r = 1:length(minRMW)
    keep = find(radRMW > minRMW(r) & radRMW <= maxRMW(r));
    length(keep);
    keep(find(keep<=40)) = [];
    for prof = 1:length(keep)
        xtmp = WSdata{keep(prof),:};
        ztmp = zdata{keep(prof),:};
        for zz = 1:length(zmax)
            keepht = find(ztmp>zmin(zz) & ztmp<=zmax(zz));
%             length(keepht)
            wspd(r,zz) = wspd(r,zz)+nansum(xtmp(keepht));
            nwspd(r,zz) = nwspd(r,zz)+length(keepht);
            clear keepht
        end
    end
    r
    clear xtmp ztmp keep
end

avwspd = wspd./nwspd; 


for r=1:24%length(minRMW)
    umaxht(r) = find(avwspd(r,:)==nanmax(avwspd(r,:)))
end

figure(2)
set(gcf,'position',[50 50 800 400])
contourf(minRMW,zmax,avwspd',100,'linecolor','none')
caxis([0 90])
shading flat,set(gca,'fontsize',14)
shg
h = colorbar;
ylabel(h,'wind speed [m s^{-1}]','fontsize',12)
ylabel('height [m]')
xlabel('R/RMW')
hold on
umaxhtmn = movmean(umaxht,9)
plot(maxRMW(1:end-1),zmax(round(umaxhtmn)),'k--')
xlim([0 4.5])

