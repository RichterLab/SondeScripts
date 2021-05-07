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
RMWmin = [0 2];
RMWmax = [2 10];
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
            if length(keep)>=5
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
        
        if length(keepfit)>5
            Ucoeffs(:,WSbin,RMWbin) = polyfit(log(meanzfit(WSbin,RMWbin,keepfit+1)+heightAdd),meanufit(WSbin,RMWbin,keepfit+1),1);
            test = fit(log(squeeze(meanzfit(WSbin,RMWbin,keepfit+1)+heightAdd)),squeeze(meanufit(WSbin,RMWbin,keepfit+1)),'poly1');
            U_ci = confint(test,0.95);
            u10(WSbin,RMWbin) = Ucoeffs(1,WSbin,RMWbin)*log(10) + Ucoeffs(2,WSbin,RMWbin);
            ustar(WSbin,RMWbin) = Ucoeffs(1,WSbin,RMWbin)*0.4;
            CD(WSbin,RMWbin) = ustar(WSbin,RMWbin)^2/u10(WSbin,RMWbin)^2;
 
            delta_u10(WSbin,RMWbin) = 2*stdufit(WSbin,RMWbin,1);
            delta_ustar(WSbin,RMWbin) = 0.5*(U_ci(2,1) - U_ci(1,1));
            delta_CD(WSbin,RMWbin) = CD(WSbin,RMWbin)*sqrt(2*(delta_ustar(WSbin,RMWbin)/abs(ustar(WSbin,RMWbin)))^2 + 2*(delta_u10(WSbin,RMWbin)/abs(u10(WSbin,RMWbin)))^2);

        else
            Ucoeffs(:,WSbin,RMWbin) = NaN;
            test = NaN;
            u10(WSbin,RMWbin) = NaN;
            ustar(WSbin,RMWbin) = NaN;
            CD(WSbin,RMWbin) = NaN;
            
            delta_u10(WSbin,RMWbin) = NaN;
            delta_ustar(WSbin,RMWbin) = NaN;
            delta_CD(WSbin,RMWbin) = NaN;
 
        end
        
    end
end

CD(find(isnan(delta_CD))) = NaN;


% Make colormap
if length(RMWmax)>5
    hex = ['#000000';'#7e1e9c';'#0165fc';'#75bbfd';'#13eac9';'#15b01a';'#9aae07';'#fac205';'#f97306';'#c65102'];
else
    hex = ['#000000';'#c65102'];
end
% 1 black, 2 ocean blue, 3 bright blue, 4 skyblue, 5 aqua
% 6 green, 7 puke green, 8 goldenrod, 9 orange, 10 dark ornage
cmap = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;

% zplot = 10:10:1e3;
% for WSbin=1:length(meanWSmin)
%     for RMWbin = 1:length(RMWmax)
%         ufitlog(WSbin,RMWbin,1:length(zplot)) = Ucoeffs(1,WSbin,RMWbin).*log(zplot) + Ucoeffs(2,WSbin,RMWbin);
%     end
% end


% Make colormap
% if length(RMWmax)>5
%     hex = ['#000000';'#7e1e9c';'#0165fc';'#75bbfd';'#13eac9';'#15b01a';'#9aae07';'#fac205';'#f97306';'#c65102'];
% else
hex = ['#000000';'#929591'];
% end
% 1 black, 2 ocean blue, 3 bright blue, 4 skyblue, 5 aqua
% 6 green, 7 puke green, 8 goldenrod, 9 orange, 10 dark ornage
cmap = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;

%% Subplot for CD vs u10
figure(1)
set(gcf,'position',[50 50 1200 700])
subplot('position',[0.55 0.56 0.43 0.4])
powellU = [27.362 33.018 40.971 50.887];
powellCD = 1e-3*[1.970 2.149 1.860 1.507];
powellust = sqrt(powellCD.*powellU.*powellU);
holtU = [19.7411, 26.6476, 33.4618,40.7699,50.9721,61.4638];
holtCD = [0.00111101,0.00161397,0.00186047,0.00212001,0.00129103,0.000744223];
holtust = sqrt(holtCD.*holtU.*holtU);
plot(powellU,powellCD,'p','markersize',10,'markerfacecolor','r','markeredgecolor','k')
hold on
plot(holtU,holtCD,'p','markersize',10,'markerfacecolor','y','markeredgecolor','k')
all = load('CD_u10_all.mat')
for wb = 1:length(meanWSmax)
    hold on
    plot(u10(wb,1),CD(wb,1),'s','color','k','markerfacecolor','k')
    plot(u10(wb,2),CD(wb,2),'s','color','r','markerfacecolor','r')
    plot(all.u10(wb),all.CD(wb),'s','color',[0.7 0.7 0.7],'markerfacecolor',[0.7 0.7 0.7])    
    errorbar(u10(wb,1),CD(wb,1),delta_CD(wb,1),'color','k')
    errorbar(u10(wb,2),CD(wb,2),delta_CD(wb,2),'color','r')
    errorbar(all.u10(wb),all.CD(wb),all.delta_CD(wb),'color',[0.7 0.7 0.7])

end
set(gca,'fontsize',14)
xlabel('{\it{U}}_{10} [m s^{-1}]')
ylim([0 0.006]),ylabel('{\it{C_D}}'),xlim([0 90])
text(3,5.5e-3,'(b)','fontsize',16)
legend('Powell et al. (2003)','Holthuijsen et al. (2012)','current study (R/RMW < 2)','current study (R/RMW > 2)','current study (all R/RMW)','box','off','location','northeast')

%% subplot of U10 vs CD with many R/RMW bins
clear all
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
RMWmin = [0 0.5 1.5 2.5 3.5 4.5 5.5 7.5];
RMWmax = [0.5 1.5 2.5 3.5 4.5 5.5 7.5 9.5];
% RMWmin = [0 2];
% RMWmax = [2 10];
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
            if length(keep)>=5
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
%         meanufit(WSbin,RMWbin,1) = NaN;
%         meanzfit(WSbin,RMWbin,1) = NaN;
        
        %% Calculate u10, CD for the near RMW region
        % Fit over 20-150 to match Vickery
        keepfit = find(~isnan(meanufit(WSbin,RMWbin,2:15)) & ~isnan(meanzfit(WSbin,RMWbin,2:15)));
        
        if length(keepfit)>5
            Ucoeffs(:,WSbin,RMWbin) = polyfit(log(meanzfit(WSbin,RMWbin,keepfit+1)+heightAdd),meanufit(WSbin,RMWbin,keepfit+1),1);
            test = fit(log(squeeze(meanzfit(WSbin,RMWbin,keepfit+1)+heightAdd)),squeeze(meanufit(WSbin,RMWbin,keepfit+1)),'poly1');
            U_ci = confint(test,0.95);
            u10(WSbin,RMWbin) = Ucoeffs(1,WSbin,RMWbin)*log(10) + Ucoeffs(2,WSbin,RMWbin);
            ustar(WSbin,RMWbin) = Ucoeffs(1,WSbin,RMWbin)*0.4;
            CD(WSbin,RMWbin) = ustar(WSbin,RMWbin)^2/u10(WSbin,RMWbin)^2;
            
            delta_u10(WSbin,RMWbin) = 2*stdufit(WSbin,RMWbin,1);
            delta_ustar(WSbin,RMWbin) = 0.5*(U_ci(2,1) - U_ci(1,1));
            delta_CD(WSbin,RMWbin) = CD(WSbin,RMWbin)*sqrt(2*(delta_ustar(WSbin,RMWbin)/abs(ustar(WSbin,RMWbin)))^2 + 2*(delta_u10(WSbin,RMWbin)/abs(u10(WSbin,RMWbin)))^2);

        else
            Ucoeffs(:,WSbin,RMWbin) = NaN;
            test = NaN;
            u10(WSbin,RMWbin) = NaN;
            ustar(WSbin,RMWbin) = NaN;
            CD(WSbin,RMWbin) = NaN;
            delta_u10(WSbin,RMWbin) = NaN;
            delta_ustar(WSbin,RMWbin) = NaN;
            delta_CD(WSbin,RMWbin) = NaN;
        end
        
    end
end

CD(find(isnan(delta_CD))) = NaN;

% Make colormap
% hex = ['#000000';'#7e1e9c';'#0165fc';'#75bbfd';'#13eac9';'#15b01a';'#9aae07';'#fac205';'#f97306';'#c65102'];
hex = ['#000000';'#7e1e9c';'#0165fc';'#13eac9';'#15b01a';'#9aae07';'#fac205';'#f97306';'#c65102'];

cmap = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;

subplot('position',[0.05 0.56 0.43 0.4])
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
    hold on
    plot(u10(wb,1),CD(wb,1),'s','color',[cmap(1,1) cmap(1,2) cmap(1,3)],'markerfacecolor',[cmap(1,1) cmap(1,2) cmap(1,3)])
    plot(u10(wb,2),CD(wb,2),'s','color',[cmap(2,1) cmap(2,2) cmap(2,3)],'markerfacecolor',[cmap(2,1) cmap(2,2) cmap(2,3)])
    plot(u10(wb,3),CD(wb,3),'s','color',[cmap(3,1) cmap(3,2) cmap(3,3)],'markerfacecolor',[cmap(3,1) cmap(3,2) cmap(3,3)])
    plot(u10(wb,4),CD(wb,4),'s','color',[cmap(4,1) cmap(4,2) cmap(4,3)],'markerfacecolor',[cmap(4,1) cmap(4,2) cmap(4,3)])
    plot(u10(wb,5),CD(wb,5),'s','color',[cmap(5,1) cmap(5,2) cmap(5,3)],'markerfacecolor',[cmap(5,1) cmap(5,2) cmap(5,3)])
    plot(u10(wb,6),CD(wb,6),'s','color',[cmap(6,1) cmap(6,2) cmap(6,3)],'markerfacecolor',[cmap(6,1) cmap(6,2) cmap(6,3)])
    plot(u10(wb,7),CD(wb,7),'s','color',[cmap(7,1) cmap(7,2) cmap(7,3)],'markerfacecolor',[cmap(7,1) cmap(7,2) cmap(7,3)])
    plot(u10(wb,8),CD(wb,8),'s','color',[cmap(8,1) cmap(8,2) cmap(8,3)],'markerfacecolor',[cmap(8,1) cmap(8,2) cmap(8,3)])
    errorbar(u10(wb,1),CD(wb,1),delta_CD(wb,1),'color',[cmap(1,1) cmap(1,2) cmap(1,3)])
    errorbar(u10(wb,2),CD(wb,2),delta_CD(wb,2),'color',[cmap(2,1) cmap(2,2) cmap(2,3)])
    errorbar(u10(wb,3),CD(wb,3),delta_CD(wb,3),'color',[cmap(3,1) cmap(3,2) cmap(3,3)])
    errorbar(u10(wb,4),CD(wb,4),delta_CD(wb,4),'color',[cmap(4,1) cmap(4,2) cmap(4,3)])
    errorbar(u10(wb,5),CD(wb,5),delta_CD(wb,5),'color',[cmap(5,1) cmap(5,2) cmap(5,3)])
    errorbar(u10(wb,6),CD(wb,6),delta_CD(wb,6),'color',[cmap(6,1) cmap(6,2) cmap(6,3)])
    errorbar(u10(wb,7),CD(wb,7),delta_CD(wb,7),'color',[cmap(7,1) cmap(7,2) cmap(7,3)])
    errorbar(u10(wb,8),CD(wb,8),delta_CD(wb,8),'color',[cmap(8,1) cmap(8,2) cmap(8,3)])
end
set(gca,'fontsize',14)
xlabel('{\it{U}}_{10} [m s^{-1}]')
ylim([0 0.006]),ylabel('{\it{C_D}}'),xlim([0 90])
text(3,5.5e-3,'(a)','fontsize',16)
hleg = legend('Powell et al. (2003)','Holthuijsen et al. (2012)','R/RMW < 0.5',...
    '0.5 < R/RMW < 1.5','1.5 < R/RMW < 2.5',...
    '2.5 < R/RMW < 3.5','3.5 < R/RMW < 4.5',...
    '4.5 < R/RMW < 5.5','5.5 < R/RMW < 7.5',...
    '7.5 < R/RMW < 10','box','off','location','northeast')
%hleg.NumColumns=2;

%% Subplot c
clearvars -except cmap
load allProfiles_3km_vmax
for n = 1:length(hurrName)
    name = hurrName{1,n};
    hurrYear(n) = str2num(name(end-3:end));
end

all_z_profiles(find(all_z_profiles==0))=NaN;
maxwind = 0;

for h = 1:length(all_U_profiles)
    ht500 = find(all_z_profiles(:,h)<=500);
    if ~isempty(ht500)
        meanWS500(h) = nanmean(all_U_profiles(ht500,h));
        maxWS500(h) = nanmax(all_U_profiles(ht500,h));
    else
        meanWS500(h) = NaN;
        maxWS500(h) = NaN;
    end
    numZ(h) = length(ht500);
    clear ht500
end

meanWSmin = [0 33 43 50 58 70]; % split by category
meanWSmax = [33 43 50 58 70 150];
% meanWSmin = [10 20 30 40 50 60 70];
% meanWSmax = [20 30 40 50 60 70 80];
RMWmin = [0 0.5 1.5 2.5 3.5 4.5 5.5 7.5];
RMWmax = [0.5 1.5 2.5 3.5 4.5 5.5 7.5 10];
RMW(find(RMW<0))=NaN;
radRMW = radius_km./RMW;
minYear = 1997;
maxYear = 2018;
heightAdd = 0;

min_fit_samples = 10;

for WSbin=1:length(meanWSmin) % 10 m/s wind speed bins
    for RMWbin = 1:length(RMWmin)
        %         clear wkeep ufit zfit meanufit meanzfit
        keep = find(meanWS500>30 & meanWS500<=40 & ...
            radRMW>RMWmin(RMWbin) & radRMW<=RMWmax(RMWbin) & ...
            stormVmax>=meanWSmin(WSbin) & stormVmax<meanWSmax(WSbin) &...
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
            if length(keep)>=5
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
        
        if length(keepfit)>5
            Ucoeffs(:,WSbin,RMWbin) = polyfit(log(meanzfit(WSbin,RMWbin,keepfit+1)+heightAdd),meanufit(WSbin,RMWbin,keepfit+1),1);
            test = fit(log(squeeze(meanzfit(WSbin,RMWbin,keepfit+1)+heightAdd)),squeeze(meanufit(WSbin,RMWbin,keepfit+1)),'poly1');
            U_ci = confint(test,0.95);
            u10(WSbin,RMWbin) = Ucoeffs(1,WSbin,RMWbin)*log(10) + Ucoeffs(2,WSbin,RMWbin);
            ustar(WSbin,RMWbin) = Ucoeffs(1,WSbin,RMWbin)*0.4;
            CD(WSbin,RMWbin) = ustar(WSbin,RMWbin)^2/u10(WSbin,RMWbin)^2;
            
            delta_u10(WSbin,RMWbin) = 2*stdufit(WSbin,RMWbin,1);
            delta_ustar(WSbin,RMWbin) = 0.5*(U_ci(2,1) - U_ci(1,1));
            delta_CD(WSbin,RMWbin) = CD(WSbin,RMWbin)*sqrt(2*(delta_ustar(WSbin,RMWbin)/abs(ustar(WSbin,RMWbin)))^2 + 2*(delta_u10(WSbin,RMWbin)/abs(u10(WSbin,RMWbin)))^2);

%             u10_low = U_ci(1,1)*log(10) + U_ci(1,2);
%             ustar_low = U_ci(1,1)*0.4;
%             CD_low(WSbin,RMWbin) = ustar_low^2/u10_low^2;
%             u10_high = U_ci(2,1)*log(10) + U_ci(2,2);
%             ustar_high = U_ci(2,1)*0.4;
%             CD_high(WSbin,RMWbin) = ustar_high^2/u10_high^2;
        else
            Ucoeffs(:,WSbin,RMWbin) = NaN;
            test = NaN;
            u10(WSbin,RMWbin) = NaN;
            ustar(WSbin,RMWbin) = NaN;
            CD(WSbin,RMWbin) = NaN;
            delta_CD(WSbin,RMWbin) = NaN;
            delta_u10(WSbin,RMWbin) = NaN;
            delta_ustar(WSbin,RMWbin) = NaN;
        end
    end
end

CD(find(isnan(delta_CD))) = NaN;

subplot('position',[0.05 0.07 0.43 0.4])
radplot = [0.25,1,2,3,4,5,6.5,8.5];
set(gca,'fontsize',12)
for cat = 1:6
    plot(radplot+0.1*cat,CD(cat,:)','s','color',[cmap(cat,1) cmap(cat,2) cmap(cat,3)],'markerfacecolor',[cmap(cat,1) cmap(cat,2) cmap(cat,3)])
    hold on
end
for cat = 1:6
    for radbin = 1:8
        errorbar(radplot(radbin)+0.1*cat,CD(cat,radbin),delta_CD(cat,radbin),'color',[cmap(cat,1) cmap(cat,2) cmap(cat,3)])
    end
end
set(gca,'fontsize',14)
xlabel('R/RMW')
ylim([0 0.006]),ylabel('{\it{C_D}}'),xlim([0 10])
legend('below category 1','category 1','category 2','category 3','category 4','category 5','box','off','location','northeast')
text(0.3,0.0055,'(c)','fontsize',16)

subplot('position',[0.55 0.07 0.43 0.4])
for cat = 1:6
    plot(radplot,numprof(cat,:)','-s','color',[cmap(cat,1) cmap(cat,2) cmap(cat,3)],'markerfacecolor',[cmap(cat,1) cmap(cat,2) cmap(cat,3)])
    hold on
end
set(gca,'fontsize',14)
xlabel('R/RMW')
ylabel('number of profiles')
legend('below category 1','category 1','category 2','category 3','category 4','category 5','box','off','location','northeast')
text(0.3,135,'(d)','fontsize',16)

