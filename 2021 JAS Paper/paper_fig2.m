close all
clear all
clc

maxwind = 0; % Set to 58 to restrict to cat 4 +, otherwise set to 0

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
% meanWSmin = [30 40 50 60 70 ];
% meanWSmax = [40 50 60 70 85 ];
% RMWmin = [0 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5];
% RMWmax = [0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5];
% RMWmin = [0 2];
% RMWmax = [2 10];
RMWmin = 0;
RMWmax = 10;
minYear = 1997;
maxYear = 2019;
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
        
        if length(keepfit)>10
            Ucoeffs(:,WSbin,RMWbin) = polyfit(log(meanzfit(WSbin,RMWbin,keepfit+1)+heightAdd),meanufit(WSbin,RMWbin,keepfit+1),1);
            test = fit(log(squeeze(meanzfit(WSbin,RMWbin,keepfit+1)+heightAdd)),squeeze(meanufit(WSbin,RMWbin,keepfit+1)),'poly1');
            U_ci = confint(test,0.95);
            u10(WSbin,RMWbin) = Ucoeffs(1,WSbin,RMWbin)*log(10) + Ucoeffs(2,WSbin,RMWbin);
            ustar(WSbin,RMWbin) = Ucoeffs(1,WSbin,RMWbin)*0.4;
            CD(WSbin,RMWbin) = ustar(WSbin,RMWbin)^2/u10(WSbin,RMWbin)^2;
            
            %compute errors:
            delta_u10(WSbin,RMWbin) = 2*stdufit(WSbin,RMWbin,1);
            delta_ustar(WSbin,RMWbin) = 0.5*(U_ci(2,1) - U_ci(1,1));
            delta_CD(WSbin,RMWbin) = CD(WSbin,RMWbin)*sqrt(2*(delta_ustar(WSbin,RMWbin)/abs(ustar(WSbin,RMWbin)))^2 + 2*(delta_u10(WSbin,RMWbin)/abs(u10(WSbin,RMWbin)))^2);
            
            
            u10_low = U_ci(1,1)*log(10) + U_ci(1,2);
            ustar_low = U_ci(1,1)*0.4;
            CD_low(WSbin,RMWbin) = ustar_low^2/u10_low^2;
            u10_high = U_ci(2,1)*log(10) + U_ci(2,2);
            ustar_high = U_ci(2,1)*0.4;
            CD_high(WSbin,RMWbin) = ustar_high^2/u10_high^2;
        else
            Ucoeffs(:,WSbin,RMWbin) = NaN;
            test = NaN;
            u10(WSbin,RMWbin) = NaN;
            ustar(WSbin,RMWbin) = NaN;
            CD(WSbin,RMWbin) = NaN;
            CD_low(WSbin,RMWbin) = NaN;
            CD_high(WSbin,RMWbin) = NaN;
        end
        
    end
end

% Make colormap
if length(RMWmax)>5
    hex = ['#000000';'#7e1e9c';'#0165fc';'#75bbfd';'#13eac9';'#15b01a';'#9aae07';'#fac205';'#f97306';'#c65102'];
else
    hex = ['#000000';'#c65102'];
end
% 1 black, 2 ocean blue, 3 bright blue, 4 skyblue, 5 aqua
% 6 green, 7 puke green, 8 goldenrod, 9 orange, 10 dark ornage
cmap = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;

%% Recreate Vickery figures 2
zplot = 10:10:1e3;
for WSbin=1:length(meanWSmin)
    for RMWbin = 1:length(RMWmax)
        ufitlog(WSbin,RMWbin,1:length(zplot)) = Ucoeffs(1,WSbin,RMWbin).*log(zplot) + Ucoeffs(2,WSbin,RMWbin);
    end
end

figure(1)
set(gcf,'position',[100 100 800 700])
subplot('position',[0.07 0.47 0.89 0.5])
for n=1:length(RMWmax) % cycle through RMW bins
    for wb = 1:length(meanWSmin)
        semilogy(squeeze(meanufit(wb,n,:)),squeeze(meanzfit(wb,n,:)),'k-','linewidth',2)
        hold on
        errorbar(squeeze(meanufit(wb,n,:)),squeeze(meanzfit(wb,n,:)),squeeze(stdufit(wb,n,:)),'horizontal','color',[0.6 0.6 0.6])
        semilogy(squeeze(ufitlog(wb,n,:)),squeeze(zplot(:)),'r','linewidth',2)
        xlim([0 100]),ylim([0 1000])
    end
    set(gca,'fontsize',14)
    set(gca,'Yscale','log')
    ylabel('height [m]','fontsize',14)
    xlabel('wind speed [m s^{-1}]','fontsize',14)
end
text(2,800,'(a)','fontsize',16)

%% Subplot for CD vs u10
subplot('position',[0.55 0.08 0.41 0.3])
powellU = [27.362 33.018 40.971 50.887];
powellCD = 1e-3*[1.970 2.149 1.860 1.507];
powellust = sqrt(powellCD.*powellU.*powellU);
holtU = [19.7411, 26.6476, 33.4618,40.7699,50.9721,61.4638];
holtCD = [0.00111101,0.00161397,0.00186047,0.00212001,0.00129103,0.000744223];
holtust = sqrt(holtCD.*holtU.*holtU);
plot(powellU,powellCD,'p','markersize',10,'markerfacecolor','r','markeredgecolor','k')
hold on
plot(holtU,holtCD,'p','markersize',10,'markerfacecolor','y','markeredgecolor','k')
for wb = 1:length(meanWSmax)
    for n = 1:length(RMWmin)
        hold on
        plot(u10(wb,n),CD(wb,n),'ks','markerfacecolor','k'); %[cmap(n,1) cmap(n,2) cmap(n,3)])
    end
end
for wb = 1:length(meanWSmax)
    for n = 1:length(RMWmin)
        hold on
        %         errorbar(u10(wb,n),CD(wb,n),abs(CD_low(wb,n)-CD(wb,n)),abs(CD_high(wb,n)-CD(wb,n)),'color','k');%[cmap(n,1) cmap(n,2) cmap(n,3)])
        errorbar(u10(wb,n),CD(wb,n),delta_CD(wb,n),'color','k');%[cmap(n,1) cmap(n,2) cmap(n,3)])
        xlim([0 100]),ylim([0 4e-3])
    end
end
set(gca,'fontsize',14)
xlabel('{\it{U}}_{10} [m s^{-1}]')
ylim([0 0.003]),ylabel('{\it{C_D}}'),xlim([0 80])
text(3,2.7e-3,'(c)','fontsize',16)

%% Subplot for u* vs u10
subplot('position',[0.07 0.08 0.41 0.3])
plot(powellU,powellust,'p','markersize',10,'markerfacecolor','r','markeredgecolor','k')
hold on
plot(holtU,holtust,'p','markersize',10,'markerfacecolor','y','markeredgecolor','k')
    plot(u10,ustar,'ks','markerfacecolor','k','markeredgecolor','k')
for wb = 1:length(meanWSmax)
        hold on
        errorbar(u10(wb,n),ustar(wb,n),delta_ustar(wb,n),'color','k');%[cmap(n,1) cmap(n,2) cmap(n,3)])
end
set(gca,'fontsize',14)
    ylabel('{\it{u_*}} [m s^{-1}]')
    xlabel('{\it{U}}_{10} [m s^{-1}]')
    ylim([0 3])

xlim([0 80])
text(3,2.7,'(b)','fontsize',16)
legend('Powell et al. (2003)','Holthuijsen et al. (2012)','current study','box','off','location','southeast')
