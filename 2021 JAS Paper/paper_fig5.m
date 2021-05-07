close all
clear all
clc

maxwind = 0; % Set to 58 to restrict to cat 4 +, otherwise set to 0
old = load('CD_u10_1997-2005.mat');

load allProfiles_3km_vmax
for n = 1:length(hurrName)
    name = hurrName{1,n};
    hurrYear(n) = str2num(name(end-3:end));
end

badData = [1253 4668 4828 5290 5296 5404 8221 6163 6287 6288 7373 9628 11434];
all_U_profiles(:,badData) = NaN; % Remove bad data, screened by Charlotte
all_z_profiles(:,badData) = NaN; % Remove bad data, screened by Charlotte
all_z_profiles(find(all_z_profiles==0))=NaN;

for h = 1:length(all_U_profiles)
    ht500 = find(all_z_profiles(:,h)<=500);
    if ~isempty(ht500)
        minWS500(h) = nanmin(all_U_profiles(ht500,h));
        meanWS500(h) = nanmean(all_U_profiles(ht500,h));
        maxWS500(h) = nanmax(all_U_profiles(ht500,h));
    else
        minWS500(h) = NaN;
        meanWS500(h) = NaN;
        maxWS500(h) = NaN;
    end
    numZ(h) = length(ht500);
    clear ht500
end

RMW(find(RMW<0))=NaN;
radRMW = radius_km./RMW;

meanWSmin = [10 20 30 40 50 60 70];
meanWSmax = [20 30 40 50 60 70 100];
RMWmin = 0;
RMWmax = 10;
minYear = 1997;
maxYear = 2005;
heightAdd = 0;

for WSbin=1:length(meanWSmin) % 10 m/s wind speed bins
    for RMWbin = 1:length(RMWmin)
        %         clear wkeep ufit zfit meanufit meanzfit
        keep = find(meanWS500>meanWSmin(WSbin) & meanWS500<=meanWSmax(WSbin) & ...
            radRMW>RMWmin(RMWbin) & radRMW<=RMWmax(RMWbin) & stormVmax>=maxwind & ...
            minZ<150 & hurrYear>=minYear & hurrYear<=maxYear);
        numprof(WSbin,RMWbin) = length(keep);
        
        ufit = nan(1,1);
        zfit = nan(1,1);
        for nearct = 1:length(keep)
            ufittmp = all_U_profiles(find(all_z_profiles(:,keep(nearct))<=1000 & all_z_profiles(:,keep(nearct))>1),keep(nearct));
            zfittmp = all_z_profiles(find(all_z_profiles(:,keep(nearct))<=1000 & all_z_profiles(:,keep(nearct))>1),keep(nearct));
            ufit = [ufit; ufittmp];
            zfit = [zfit; zfittmp];
        end
        
        for ht = 1:100 % 10-m height bins
            %% Do the near ones
            keep = find(zfit>(ht-1)*10 &  zfit<=ht*10);
            numpts(WSbin,RMWbin,ht) = length(keep);
            if length(keep)>=10
                meanufit(WSbin,RMWbin,ht) = nanmean(ufit(keep));
                meanzfit(WSbin,RMWbin,ht) = nanmean(zfit(keep));
                stdufit(WSbin,RMWbin,ht) = nanstd(ufit(keep));
            else
                meanufit(WSbin,RMWbin,ht) = NaN;
                meanzfit(WSbin,RMWbin,ht) = NaN;
                stdufit(WSbin,RMWbin,ht) = NaN;
            end
        end
        
        % Exclude the lowest 20 m because the log layer may not be valid there
        meanufit(WSbin,RMWbin,1) = NaN;
        meanzfit(WSbin,RMWbin,1) = NaN;
    end
end

%% Recreate Vickery figures 2
zplot = 10:10:1e3;
for WSbin=1:length(meanWSmin)
    for RMWbin = 1:length(RMWmax)
        ufitlog(WSbin,RMWbin,1:length(zplot)) = old.Ucoeffs(1,WSbin,RMWbin).*log(zplot) + old.Ucoeffs(2,WSbin,RMWbin);
    end
end

figure(1)
set(gcf,'position',[100 100 1200 720])
subplot('position',[0.05 0.51 0.43 0.46])
for n=1:length(RMWmax) % cycle through RMW bins
    for wb = 1:length(meanWSmin)
        semilogy(squeeze(meanufit(wb,n,:)),squeeze(meanzfit(wb,n,:)),'k-','linewidth',2)
        hold on
        errorbar(squeeze(meanufit(wb,n,:)),squeeze(meanzfit(wb,n,:)),squeeze(stdufit(wb,n,:)),'horizontal','color',[0.6 0.6 0.6])
        semilogy(squeeze(ufitlog(wb,n,:)),squeeze(zplot(:)),'r','linewidth',2)
    end
    set(gca,'fontsize',14)
    set(gca,'Yscale','log')
    ylabel('height [m]','fontsize',14)
    xlabel('wind speed [m s^{-1}]','fontsize',14)
end
% set(gcf, 'DefaultTextBackgroundColor', [1,1,1])
for wb = 1:length(meanWSmin)
    text(ufitlog(wb,n,50-wb*4)+0.5,zplot(50-wb*4),num2str(old.numprof(wb)),'fontsize',16)
end

text(2,800,'(a)','fontsize',16)
xlim([0 90]),ylim([0 1000])

clear all

maxwind = 0;

load allProfiles_3km_vmax
new = load('CD_u10_2006-2018.mat');

for n = 1:length(hurrName)
    name = hurrName{1,n};
    hurrYear(n) = str2num(name(end-3:end));
end

badData = [1253 4668 4828 5290 5296 5404 8221 6163 6287 6288 7373 9628 11434];
all_U_profiles(:,badData) = NaN; % Remove bad data, screened by Charlotte
all_z_profiles(:,badData) = NaN; % Remove bad data, screened by Charlotte
all_z_profiles(find(all_z_profiles==0))=NaN;

for h = 1:length(all_U_profiles)
    ht500 = find(all_z_profiles(:,h)<=500);
    if ~isempty(ht500)
        minWS500(h) = nanmin(all_U_profiles(ht500,h));
        meanWS500(h) = nanmean(all_U_profiles(ht500,h));
        maxWS500(h) = nanmax(all_U_profiles(ht500,h));
    else
        minWS500(h) = NaN;
        meanWS500(h) = NaN;
        maxWS500(h) = NaN;
    end
    numZ(h) = length(ht500);
    clear ht500
end

RMW(find(RMW<0))=NaN;
radRMW = radius_km./RMW;

meanWSmin = [10 20 30 40 50 60 70];
meanWSmax = [20 30 40 50 60 70 100];
RMWmin = 0;
RMWmax = 10;
minYear = 2006;
maxYear = 2018;
heightAdd = 0;

for WSbin=1:length(meanWSmin) % 10 m/s wind speed bins
    for RMWbin = 1:length(RMWmin)
        %         clear wkeep ufit zfit meanufit meanzfit
        keep = find(meanWS500>meanWSmin(WSbin) & meanWS500<=meanWSmax(WSbin) & ...
            radRMW>RMWmin(RMWbin) & radRMW<=RMWmax(RMWbin) & stormVmax>=maxwind & ...
            minZ<150 & hurrYear>=minYear & hurrYear<=maxYear);
        numprof(WSbin,RMWbin) = length(keep);
        
        ufit = nan(1,1);
        zfit = nan(1,1);
        for nearct = 1:length(keep)
            ufittmp = all_U_profiles(find(all_z_profiles(:,keep(nearct))<=1000 & all_z_profiles(:,keep(nearct))>1),keep(nearct));
            zfittmp = all_z_profiles(find(all_z_profiles(:,keep(nearct))<=1000 & all_z_profiles(:,keep(nearct))>1),keep(nearct));
            ufit = [ufit; ufittmp];
            zfit = [zfit; zfittmp];
        end
        
        for ht = 1:100 % 10-m height bins
            %% Do the near ones
            keep = find(zfit>(ht-1)*10 &  zfit<=ht*10);
            numpts(WSbin,RMWbin,ht) = length(keep);
            if length(keep)>=10
                meanufit(WSbin,RMWbin,ht) = nanmean(ufit(keep));
                meanzfit(WSbin,RMWbin,ht) = nanmean(zfit(keep));
                stdufit(WSbin,RMWbin,ht) = nanstd(ufit(keep));
            else
                meanufit(WSbin,RMWbin,ht) = NaN;
                meanzfit(WSbin,RMWbin,ht) = NaN;
                stdufit(WSbin,RMWbin,ht) = NaN;
            end
        end
        
        % Exclude the lowest 20 m because the log layer may not be valid there
        meanufit(WSbin,RMWbin,1) = NaN;
        meanzfit(WSbin,RMWbin,1) = NaN;
        
        %% Calculate u10, CD for the near RMW region
        % Fit over 20-150 to match Vickery
        keepfit = find(~isnan(meanufit(WSbin,RMWbin,2:15)) & ~isnan(meanzfit(WSbin,RMWbin,2:15)));
        
    end
end


%% Recreate Vickery figures 2
zplot = 10:10:1e3;
for WSbin=1:length(meanWSmin)
    for RMWbin = 1:length(RMWmax)
        ufitlog(WSbin,RMWbin,1:length(zplot)) = new.Ucoeffs(1,WSbin,RMWbin).*log(zplot) + new.Ucoeffs(2,WSbin,RMWbin);
    end
end

subplot('position',[0.54 0.51 0.43 0.46])
for n=1:length(RMWmax) % cycle through RMW bins
    for wb = 1:length(meanWSmin)
        semilogy(squeeze(meanufit(wb,n,:)),squeeze(meanzfit(wb,n,:)),'k-','linewidth',2)
        hold on
        errorbar(squeeze(meanufit(wb,n,:)),squeeze(meanzfit(wb,n,:)),squeeze(stdufit(wb,n,:)),'horizontal','color',[0.6 0.6 0.6])
        semilogy(squeeze(ufitlog(wb,n,:)),squeeze(zplot(:)),'r','linewidth',2)
        xlim([0 90]),ylim([0 1000])
    end
    set(gca,'fontsize',14)
    set(gca,'Yscale','log')
    ylabel('height [m]','fontsize',14)
    xlabel('wind speed [m s^{-1}]','fontsize',14)
end
for wb = 1:length(meanWSmin)
    text(ufitlog(wb,n,50-wb*4)+0.5,zplot(50-wb*4),num2str(new.numprof(wb)),'fontsize',16)
end
text(2,800,'(b)','fontsize',16)

%% Subplot for CD vs u10
old = load('CD_u10_1997-2005.mat');
new = load('CD_u10_2006-2018.mat');


subplot('position',[0.05 0.09 0.92 0.33])
powellU = [27.362 33.018 40.971 50.887];
powellCD = 1e-3*[1.970 2.149 1.860 1.507];
powellust = sqrt(powellCD.*powellU.*powellU);
holtU = [19.7411, 26.6476, 33.4618,40.7699,50.9721,61.4638];
holtCD = [0.00111101,0.00161397,0.00186047,0.00212001,0.00129103,0.000744223];
holtust = sqrt(holtCD.*holtU.*holtU);
plot(powellU,powellCD,'p','markersize',10,'markerfacecolor','r','markeredgecolor','k')
hold on
plot(holtU,holtCD,'p','markersize',10,'markerfacecolor','y','markeredgecolor','k')
plot(old.u10(1),old.CD(1),'ks','markerfacecolor','k'); %[cmap(n,1) cmap(n,2) cmap(n,3)])
plot(new.u10(1),new.CD(1),'s','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor',[0.7 0.7 0.7]); %[cmap(n,1) cmap(n,2) cmap(n,3)])


for wb = 1:length(meanWSmax)
    for n = 1:length(RMWmin)
        hold on
        plot(old.u10(wb,n),old.CD(wb,n),'ks','markerfacecolor','k'); %[cmap(n,1) cmap(n,2) cmap(n,3)])
        errorbar(old.u10(wb,n),old.CD(wb,n),old.delta_CD(wb,n),'color','k');%[cmap(n,1) cmap(n,2) cmap(n,3)])
        plot(new.u10(wb,n),new.CD(wb,n),'s','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor',[0.7 0.7 0.7]); %[cmap(n,1) cmap(n,2) cmap(n,3)])
        errorbar(new.u10(wb,n),new.CD(wb,n),new.delta_CD(wb,n),'color',[0.7 0.7 0.7]);%[cmap(n,1) cmap(n,2) cmap(n,3)])
    end
end
set(gca,'fontsize',14)
text(1.25,4.5e-3,'(c)','fontsize',16)
xlabel('{\it{U}}_{10} [m s^{-1}]')
ylim([0 0.005]),ylabel('{\it{C_D}}'),xlim([0 80])

load pick3
CDvals = CDsv;
maxCD = nanmax(CDvals');
minCD = nanmin(CDvals');

ustarvals = ustsv;
maxustar = nanmax(ustarvals');
minustar = nanmin(ustarvals');

u10vals = u10sv;
u10 = nanmedian(u10vals');
% fill([u10 fliplr(u10)],[minCD fliplr(maxCD)],[0.4 0.9 0.3],'facealpha',0.5,'edgecolor','none')
for wb = 1:length(meanWSmax)
    for n = 1:length(RMWmin)
        hold on
        plot(old.u10(wb,n),old.CD(wb,n),'ks','markerfacecolor','k'); %[cmap(n,1) cmap(n,2) cmap(n,3)])
        errorbar(old.u10(wb,n),old.CD(wb,n),old.delta_CD(wb,n),'color','k');%[cmap(n,1) cmap(n,2) cmap(n,3)])
        plot(new.u10(wb,n),new.CD(wb,n),'s','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor',[0.7 0.7 0.7]); %[cmap(n,1) cmap(n,2) cmap(n,3)])
        errorbar(new.u10(wb,n),new.CD(wb,n),new.delta_CD(wb,n),'color',[0.7 0.7 0.7]);%[cmap(n,1) cmap(n,2) cmap(n,3)])
    end
end
plot(powellU,powellCD,'p','markersize',10,'markerfacecolor','r','markeredgecolor','k')
plot(holtU,holtCD,'p','markersize',10,'markerfacecolor','y','markeredgecolor','k')

legend('Powell et al. (2003)','Holthuijsen et al. (2012)','current study (1997 - 2005)','current study (2006 - 2018)','box','off','location','northeast')
