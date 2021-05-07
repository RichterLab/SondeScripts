close all
clear all
clc

%% Make figure 9 for the paper
x = squeeze(ncread('azim_avg_new.nc','xh'));
z = squeeze(ncread('azim_avg_new.nc','zh'));
u = squeeze(ncread('azim_avg_new.nc','uinterp'));
v = squeeze(ncread('azim_avg_new.nc','vinterp'));
ust = squeeze(ncread('azim_avg_new.nc','ust'));
s10 = squeeze(ncread('azim_avg_new.nc','s10'));
RMW = 12e3;
radRMW = x./RMW;
radstart = 5;

WS = sqrt(u.^2 + v.^2);
[Z X] = meshgrid(z,x);

[nr nz] = size(u);

cmap = jet(10);
hex = ['#929591';'#7e1e9c';'#0165fc';'#75bbfd';'#13eac9';'#15b01a';'#9aae07';'#fac205';'#f97306';'#c65102'];
cmap = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;


figure(1)
set(gcf,'position',[100 100 800 700])
hold all
cmapct = 1;
for ir = radstart:100:find(x<30000,1,'last')
    subplot('position',[0.07 0.61 0.89 0.36])
    hold on
    plot(WS(ir,1:64),Z(ir,1:64),'-s','color',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)],'markerfacecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)])
    
    Ucoeffs = polyfit(log(z(2:4)),WS(ir,2:4),1) % fit from 23 to 100 m
    uline = polyval(Ucoeffs,log(z(2:4)));
    %     plot(uline,z(2:4),'k')
    set(gca,'yscale','log')
    set(gca,'fontsize',14)
    ylabel('height [m]')
    xlabel('wind speed [m s^{-1}]')
    xlim([0 120])
    
    u10(cmapct) = Ucoeffs(1)*log(10) + Ucoeffs(2);
    ustar(cmapct) = Ucoeffs(1)*0.4;
    CD(cmapct) = ustar(cmapct).^2/u10(cmapct).^2;
    clear Ucoeffs ufit zfit R
    
    subplot('position',[0.07 0.335 0.42 0.2])
    plot(u10(cmapct),CD(cmapct),'s','markerfacecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)],'markeredgecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)])
    hold on
    ylabel('\it{C_D}','fontsize',14)
    xlabel('{\it{U}}_{10} [m s^{-1}]')
    set(gca,'fontsize',14)
    
    subplot('position',[0.07 0.06 0.42 0.2])
    plot(radRMW(ir),CD(cmapct),'s','markerfacecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)],'markeredgecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)])
    hold on
    ylabel('\it{C_D}','fontsize',14)
    xlabel('R/RMW')
    set(gca,'fontsize',14)
    
    subplot('position',[0.55 0.335 0.42 0.2])
    plot(u10(cmapct),ustar(cmapct),'s','markerfacecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)],'markeredgecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)])
    hold on
    set(gca,'fontsize',14)
    ylabel('{\it{u_*}} [m s^{-1}]','fontsize',14)
    xlabel('{\it{U}}_{10} [m s^{-1}]')
    
    subplot('position',[0.55 0.06 0.42 0.2])
    plot(radRMW(ir),ustar(cmapct),'s','markerfacecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)],'markeredgecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)])
    hold on
    ylabel('{\it{u_*}} [m s^{-1}]','fontsize',14)
    xlabel('R/RMW')
    set(gca,'fontsize',14)    
    
    radius{cmapct} = ['R/RMW = ' num2str(1e-3*x(ir)/12,'%.2f')];
    
    cmapct = cmapct+1;
end

extra = load('fig10extra.mat');

subplot('position',[0.07 0.61 0.89 0.36])
legend(radius,'box','off')

subplot('position',[0.07 0.61 0.89 0.36])
text(3,0.7e3,'(a)','fontsize',16)

subplot('position',[0.07 0.335 0.42 0.2])
text(3,2.6e-3,'(b)','fontsize',16)
ufix = [0 25 70];
cdfix = [1e-3 2.4e-3 2.4e-3];
ustfix = sqrt(cdfix.*ufix.*ufix)
plot(ufix,cdfix,'k--')
plot(extra.u10,extra.CD,'o','markerfacecolor',[0.8 0.8 0.8],'markeredgecolor',[0.8 0.8 0.8])
plot(u10,CD,'-','color',[0.5 0.5 0.5])


subplot('position',[0.07 0.06 0.42 0.2])
text(0.1,2.6e-3,'(d)','fontsize',16)
plot(extra.radRMW(1:10),extra.CD,'o','markerfacecolor',[0.8 0.8 0.8],'markeredgecolor',[0.8 0.8 0.8])
plot(radRMW(1:100:901),CD,'-','color',[0.5 0.5 0.5])

subplot('position',[0.55 0.335 0.42 0.2])
text(3,3.45,'(c)','fontsize',16)
plot(ufix,ustfix,'k--')
plot(extra.u10,extra.ustar,'o','markerfacecolor',[0.8 0.8 0.8],'markeredgecolor',[0.8 0.8 0.8])
plot(u10,ustar,'-','color',[0.5 0.5 0.5])


subplot('position',[0.55 0.06 0.42 0.2])
text(0.1,3.45,'(e)','fontsize',16)
plot(extra.radRMW(1:10),extra.ustar,'o','markerfacecolor',[0.8 0.8 0.8],'markeredgecolor',[0.8 0.8 0.8])
plot(radRMW(1:100:901),ustar,'-','color',[0.5 0.5 0.5])

cmapct = 1;
for ir = radstart:100:find(x<30000,1,'last')
    
%     u10(cmapct) = Ucoeffs(1)*log(10) + Ucoeffs(2);
%     ustar(cmapct) = Ucoeffs(1)*0.4;
%     CD(cmapct) = ustar(cmapct).^2/u10(cmapct).^2;
%     clear Ucoeffs ufit zfit R
    
    subplot('position',[0.07 0.06 0.42 0.2])
    plot(radRMW(ir),CD(cmapct),'s','markerfacecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)],'markeredgecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)])
    hold on
    ylabel('\it{C_D}','fontsize',14)
    xlabel('R/RMW')
    set(gca,'fontsize',14)
       
    subplot('position',[0.55 0.06 0.42 0.2])
    plot(radRMW(ir),ustar(cmapct),'s','markerfacecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)],'markeredgecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)])
    hold on
    ylabel('{\it{u_*}} [m s^{-1}]','fontsize',14)
    xlabel('R/RMW')
    set(gca,'fontsize',14)    
    
%     radius{cmapct} = ['R/RMW = ' num2str(1e-3*x(ir)/12,'%.2f')];
    
    cmapct = cmapct+1;
end


hold all
cmapct = 1;
for ir = radstart:100:find(x<30000,1,'last')
    subplot('position',[0.07 0.61 0.89 0.36])
    hold on
    Ucoeffs = polyfit(log(z(2:4)),WS(ir,2:4),1) % fit from 23 to 100 m
    uline = polyval(Ucoeffs,log(z(2:4)));
    plot(uline,z(2:4),'k','linewidth',2,'HandleVisibility','off')
    
    u10(cmapct) = Ucoeffs(1)*log(10) + Ucoeffs(2);
    ustar(cmapct) = Ucoeffs(1)*0.4;
    CD(cmapct) = ustar(cmapct).^2/u10(cmapct).^2;
    clear Ucoeffs ufit zfit R
    
        subplot('position',[0.07 0.335 0.42 0.2])
    plot(u10(cmapct),CD(cmapct),'s','markerfacecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)],'markeredgecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)])
    hold on
    ylabel('\it{C_D}','fontsize',14)
    xlabel('{\it{U}}_{10} [m s^{-1}]')
    set(gca,'fontsize',14)
    
        subplot('position',[0.55 0.335 0.42 0.2])
    plot(u10(cmapct),ustar(cmapct),'s','markerfacecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)],'markeredgecolor',[cmap(cmapct,1) cmap(cmapct,2) cmap(cmapct,3)])
    hold on
    set(gca,'fontsize',14)
    ylabel('{\it{u_*}} [m s^{-1}]','fontsize',14)
    xlabel('{\it{U}}_{10} [m s^{-1}]')
    
    cmapct = cmapct+1;
end


