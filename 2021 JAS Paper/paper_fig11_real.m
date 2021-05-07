close all
clear all
clc

load windProfiles_allSondes.mat

%% Try dynamic fitting over a range of heights to see where the log layer
% assumption breaks down

% RMWmin = 0.75;
% RMWmax = 1.25;

for ww = 1:size(meanufit,1) % wind speed bin
    if ~isnan(meanufit(ww,3))
        % start fitting over the lowest 30 m and go from there
        % height bin is 5 m
        for maxht = 10:30
            ufit = squeeze(meanufit(ww,2:maxht));
            zfit = squeeze(meanzfit(ww,2:maxht));
            
            R = corrcoef(ufit,log(zfit));
            Rsq(ww,maxht-9) = R(1,2).^2;
            
            Ucoeffs = polyfit(log(zfit),ufit,1);
            u10(ww,maxht-9) = Ucoeffs(1)*log(10) + Ucoeffs(2);
            ustar(ww,maxht-9) = Ucoeffs(1)*0.4;
            CD(ww,maxht-9) = ustar(ww,maxht-9).^2/u10(ww,maxht-9).^2;
            clear Ucoeffs ufit zfit
        end
    end
end

Rsq(find(Rsq==0))=NaN;
CD(find(CD==0))=NaN;
u10(find(u10==0))=NaN;
ustar(find(ustar==0))=NaN;

ufix = [0 25 90];
cdfix = [1e-3 2.4e-3 2.4e-3];

%% Test how height varies CD in depth for one bin
% R/RMW = 1.5-1.75 (7), wind speed = 65-70 m/s (12)

% start fitting over the lowest 30 m and go from there
% height bin is 5 m
windbin = 6;
figure(11)
set(gcf,'position',[100 100 1000 720])
subplot('position',[0.06 0.08 0.37 0.88])
semilogy(squeeze(meanufit(windbin,1,2:end)),squeeze(meanzfit(windbin,1,2:end)),'k','linewidth',2)
set(gca,'fontsize',12)
% title(['R/RMW = 0.5-1.5, WS_{500} = ' num2str(meanWSmin(windbin)) '-' num2str(meanWSmax(windbin)) ' m s^{-1}'])
hold on
cscl = jet(21);
for maxht = 10:30
    ufit = squeeze(meanufit(windbin,1,2:maxht));
    zfit = squeeze(meanzfit(windbin,1,2:maxht));
    
    R = corrcoef(ufit,log(zfit));
    Rsq_test(maxht-9) = R(1,2).^2;
    
    Ucoeffs = polyfit(log(zfit),ufit,1);
    uline = polyval(Ucoeffs,log(zfit));
    subplot('position',[0.06 0.08 0.37 0.88])
    semilogy(uline,zfit,'color',[cscl(maxht-9,1) cscl(maxht-9,2) cscl(maxht-9,3)])
    ylabel('{\it{z}} [m]'),ylim([10 500]),
    %     xlim([60 68])
    xlabel('wind speed [m s^{-1}]')
    set(gca,'fontsize',14)
    set(gca,'ytick',[10 20 50 100 200 500])
    
    u10_test(maxht-9) = Ucoeffs(1)*log(10) + Ucoeffs(2);
    ustar_test(maxht-9) = Ucoeffs(1)*0.4;
    CD_test(maxht-9) = ustar_test(maxht-9).^2/u10_test(maxht-9).^2;
    
    subplot('position',[0.5 0.71 0.46 0.24])
    semilogy(Rsq_test(maxht-9),nanmax(zfit),'s','linewidth',2,'color',[cscl(maxht-9,1) cscl(maxht-9,2) cscl(maxht-9,3)],'markerfacecolor',[cscl(maxht-9,1) cscl(maxht-9,2) cscl(maxht-9,3)])
    hold on
    ylim([100 300]),set(gca,'fontsize',12)
    xlabel('{\it{R}}^2','fontsize',14)
    ylabel('Maximum fitting height [m]')
    set(gca,'fontsize',14)
    
    subplot('position',[0.5 0.395 0.46 0.24])
    semilogy(ustar_test(maxht-9),nanmax(zfit),'s','linewidth',2,'color',[cscl(maxht-9,1) cscl(maxht-9,2) cscl(maxht-9,3)],'markerfacecolor',[cscl(maxht-9,1) cscl(maxht-9,2) cscl(maxht-9,3)])
    hold on
    ylim([100 300]),set(gca,'fontsize',12),xlabel('{\it{u_*}} [m s^{-1}]')
    ylabel('Maximum fitting height [m]')
    set(gca,'fontsize',14)
    
    subplot('position',[0.5 0.08 0.46 0.24])
    semilogy(CD_test(maxht-9),nanmax(zfit),'s','linewidth',2,'color',[cscl(maxht-9,1) cscl(maxht-9,2) cscl(maxht-9,3)],'markerfacecolor',[cscl(maxht-9,1) cscl(maxht-9,2) cscl(maxht-9,3)])
    hold on
    ylim([100 300]),set(gca,'fontsize',12),xlabel('{\it{C_D}}')
    ylabel('Maximum fitting height [m]')
    set(gca,'fontsize',14)
    
    clear Ucoeffs ufit zfit
end

subplot('position',[0.06 0.08 0.37 0.88])
text(68,430,'(a)','fontsize',16)
semilogy(squeeze(meanufit(windbin,1,2:end)),squeeze(meanzfit(windbin,1,2:end)),'k','linewidth',2)

subplot('position',[0.5 0.71 0.46 0.24])
text(0.996,270,'(b)','fontsize',16)
xlim([0.95 1])

subplot('position',[0.5 0.395 0.46 0.24])
text(2.17,270,'(c)','fontsize',16)
xlim([1.8 2.2])

subplot('position',[0.5 0.08 0.46 0.24])
text(1.77e-3,270,'(d)','fontsize',16)
xlim([1.4e-3 1.8e-3])