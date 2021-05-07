close all
clear all
clc

% meds = load('fit_z10_z150_median.mat');
load pick3
z10z150 = load('fit_z10_z150.mat');


CDvals = CDsv;
maxCD = nanmax(CDvals');
minCD = nanmin(CDvals');

ustarvals = ustsv;
maxustar = nanmax(ustarvals');
minustar = nanmin(ustarvals');

u10vals = u10sv;
u10 = nanmedian(u10vals');

figure(1)
set(gcf,'position',[100 100 800 300])
subplot('position',[0.07 0.18 0.41 0.76])
fill([u10 fliplr(u10)],[minCD fliplr(maxCD)],[0.4 0.9 0.3],'facealpha',0.5,'edgecolor','none')
hold on
plot(z10z150.u10,z10z150.CD,'k-s','linewidth',2,'markerfacecolor','k')
set(gca,'fontsize',14)
ylabel('{\it{C_D}}')
xlabel('{\it{U}}_{10} [m s^{-1}]')
ylim([0 3e-3]),xlim([0 60])
text(3,2.7e-3,'(a)','fontsize',14)
    
subplot('position',[0.55 0.18 0.41 0.76])
fill([z10z150.u10;flipud(z10z150.u10)],[minustar fliplr(maxustar)],[0.4 0.9 0.3],'facealpha',0.5,'edgecolor','none')
hold on
plot(z10z150.u10,z10z150.ustar,'k-s','linewidth',2,'markerfacecolor','k')
set(gca,'fontsize',14)
ylabel('{\it{u_*}} [m s^{-1}]')
xlabel('{\it{U}}_{10} [m s^{-1}]')
ylim([0 3]),xlim([0 60])
text(3,2.7,'(b)','fontsize',14)


% %% Subplot for CD vs u10
% subplot('position',[0.07 0.08 0.41 0.3])
% powellU = [27.362 33.018 40.971 50.887];
% powellCD = 1e-3*[1.970 2.149 1.860 1.507];
% powellust = sqrt(powellCD.*powellU.*powellU);
% holtU = [19.7411, 26.6476, 33.4618,40.7699,50.9721,61.4638];
% holtCD = [0.00111101,0.00161397,0.00186047,0.00212001,0.00129103,0.000744223];
% holtust = sqrt(holtCD.*holtU.*holtU);
% plot(powellU,powellCD,'p','markersize',10,'markerfacecolor','r','markeredgecolor','k')
% hold on
% plot(holtU,holtCD,'p','markersize',10,'markerfacecolor','y','markeredgecolor','k')
% for wb = 1:length(meanWSmax)
%     for n = 1:length(RMWmin)
%         hold on
%         plot(u10(wb,n),CD(wb,n),'ks','markerfacecolor','k'); %[cmap(n,1) cmap(n,2) cmap(n,3)])
%     end
% end
% for wb = 1:length(meanWSmax)
%     for n = 1:length(RMWmin)
%         hold on
%         errorbar(u10(wb,n),CD(wb,n),abs(CD_low(wb,n)-CD(wb,n)),abs(CD_high(wb,n)-CD(wb,n)),'color','k');%[cmap(n,1) cmap(n,2) cmap(n,3)])
%         xlim([0 100]),ylim([0 4e-3])
%     end
% end
% set(gca,'fontsize',14)
% xlabel('{\it{U}}_{10} [m s^{-1}]')
% ylim([0 0.003]),ylabel('{\it{C_D}}'),xlim([0 80])
% text(3,2.7e-3,'(b)','fontsize',16)
% 
% %% Subplot for u* vs u10
% subplot('position',[0.55 0.08 0.41 0.3])
% plot(powellU,powellust,'p','markersize',10,'markerfacecolor','r','markeredgecolor','k')
% hold on
% plot(holtU,holtust,'p','markersize',10,'markerfacecolor','y','markeredgecolor','k')
% for n=1:length(RMWmax) % cycle through RMW bins
%     plot(u10,ustar,'ks','markerfacecolor','k','markeredgecolor','k')
%     set(gca,'fontsize',14)
%     ylabel('{\it{u_*}} [m s^{-1}]')
%     xlabel('{\it{U}}_{10} [m s^{-1}]')
%     ylim([0 3])
% end
% xlim([0 80])
% text(3,2.7,'(c)','fontsize',16)
% legend('Powell et al. (2003)','Holthuijsen et al. (2012)','current study','box','off','location','southeast')
% 
