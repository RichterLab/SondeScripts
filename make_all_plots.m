%% 1. Read the output data(from get_velocity_h_r_w.m) for different height (100:100:500) and radius (R/RMW)
% figure
% count=0;
for Fitting_height=100:100:500
clc;
clearvars -except count u10* CD* ustar* Fitting_height
%Split into 10m/s-wide bins from 0 to 80 m/s (each with its own mean profile)
%Split vertical into 5m-wide bins
% count=count+1

for num_z_define=398;%20:10:50

Fitting_name=strcat('allStorms_constants_rad_MBL',num2str(Fitting_height),'.mat')
load(Fitting_name);

min_samples = 8;   %Minimum number of samples to be included
min_fit_samples = 8;
%build height bin
zplot = min_height+0.5*height_interval:height_interval:max_height-0.5*height_interval;

max_wind = 80;                             
wind_interval = 10;
num_wind_bins = max_wind/wind_interval;

start_fitting=2;%20m
end_fitting=Fitting_height/height_interval-start_fitting;%100m

%%plot
%All includes 2425 profiles
Fitting_name=strcat('allStorms_all_profile_data_rad_MBL',num2str(Fitting_height),'.mat')
load(Fitting_name);
num_z=num_z_define;
%Now fit the log region in these mean profiles:
Ucoeffs = zeros(2,num_rad_bins,num_wind_bins);
ustar   = zeros(num_rad_bins,num_wind_bins);
u10     = zeros(num_rad_bins,num_wind_bins);
CD      = zeros(num_rad_bins,num_wind_bins);

mean_U_profiles = zeros(num_z,num_rad_bins,num_wind_bins);
std_U_profiles  = zeros(num_z,num_rad_bins,num_wind_bins);
numvecU         = zeros(num_z,num_rad_bins,num_wind_bins);

for j=1:num_wind_bins % wind radius
    for i=1:num_rad_bins % bin radius
        for k=1:num_z % bin height
            mean_U_profiles(k,i,j) = mean(all_U_profiles{k,i,j}(:));
            std_U_profiles(k,i,j) = std(all_U_profiles{k,i,j}(:));
            numvecU(k,i,j) = length(all_U_profiles{k,i,j}(:));
        end %bin height
    end % bin radius
end %wind velocity
%
count_fitting=0;
for num_fitting=13:5:63
    count_fitting=count_fitting+1;
    for j=1:num_wind_bins % wind radius
        for i=1:num_rad_bins % bin radius
            tmp = mean_U_profiles(3:num_fitting,i,j);
            numtmp = numvecU(3:num_fitting,i,j);
            zfit = zplot(~isnan(tmp)&numtmp>min_samples); %what are we removing here?
            ufit = tmp(~isnan(tmp)&numtmp>min_samples);
            if (length(ufit) > min_fit_samples)        
                [curve2,gof2] = fit(log(zfit)',ufit,'poly1');
                err_fit(i,j,count_fitting)=gof2.rmse;
            end
        end
    end
end

for j=1:num_wind_bins % wind radius
    for i=1:num_rad_bins % bin radius
    num_fitting=12;
    tmp = mean_U_profiles(start_fitting:end_fitting,i,j);
    numtmp = numvecU(start_fitting:end_fitting,i,j);
    zfit = zplot(~isnan(tmp)&numtmp>min_samples); %what are we removing here?
    ufit = tmp(~isnan(tmp)&numtmp>min_samples);
    if (length(ufit) > min_fit_samples)
        Ucoeffs(:,i,j) = polyfit(log(zfit),ufit',1); % why log(x) fit?
        test = fit(log(zfit)',ufit,'poly1');
        U_ci = confint(test,0.95);
    else
        Ucoeffs(:,i,j) = ones(2,1)*NaN;
        U_ci = ones(2,2)*NaN;
    end
    
    %Compute the mean velocity, temperature near the 10-meter height:
    kz10 = floor(10-min_height)/height_interval + 1;  %only if min_height < 10
    z10 = zplot(kz10);
    u10(i,j) = Ucoeffs(1,i,j)*log(10) + Ucoeffs(2,i,j);
    
    ustar(i,j) = Ucoeffs(1,i,j)*0.4;
    CD(i,j) = ustar(i,j)^2/u10(i,j)^2;
    
    %compute errors:
    delta_u10(i,j) = 2*std_U_profiles(1,i,j);
    delta_ustar(i,j) = 0.5*(U_ci(2,1) - U_ci(1,1));
    delta_CD(i,j) = CD(i,j)*sqrt(2*(delta_ustar(i,j)/abs(ustar(i,j)))^2 + 2*(delta_u10(i,j)/abs(u10(i,j)))^2);
    
    end
end

% % subplot(2,2,count)
% % for i=1:1:4
% %     hold on
% %     p(i)=plot(u10(i,:),CD(i,:),'-*');
% % end
% % set(p(1),'linewidth',2.0,'color','k');
% % set(p(2),'linewidth',2.0,'color','b');
% % set(p(3),'linewidth',2.0,'color','r');
% % set(p(4),'linewidth',2.0,'color','g');
% % % set(p(5),'linewidth',2.0,'color','c');
% % % set(p(6),'linewidth',2.0,'color','y');
% % % 
% % % set(p(1+6),'linewidth',2.0,'color','k');
% % % set(p(2+6),'linewidth',2.0,'color','b');
% % % set(p(3+6),'linewidth',2.0,'color','r');
% % % set(p(4+6),'linewidth',2.0,'color','g');
% % % set(p(5+6),'linewidth',2.0,'color','c');
% % % set(p(6+6),'linewidth',2.0,'color','y');
% % 
% % %
% % u10plt = 0:2:70;
% % LP_CD = (0.49 + 0.065.*u10plt)./1000;
% % LP_ustar = sqrt(LP_CD.*u10plt.^2);
% % Andreas_ustar = 0.0583*u10plt - 0.243;
% % Andreas_CD = Andreas_ustar.^2./u10plt.^2;
% % powellCDsquares = load('./Published_data/powellCDsquares.dat');
% % Powellustar = load('./Published_data/Powell_ustar.dat');
% % load('./Published_data/ZhangData.mat');
% % HolthuijsenCD = load('./Published_data/HolthuijsenCD.dat');
% % HolthuijsenCD_errortop = load('./Published_data/Holthuijsen_CD_errortop.dat');
% % HolthuijsenCD_error = (HolthuijsenCD_errortop(:,2)-HolthuijsenCD(:,2));
% % 
% % hold on
% % % errorbar(u10(:),CD(:),delta_CD(:),'sb','markerfacecolor','b')
% % p(13)= plot(u10plt(:),LP_CD(:),'--k','linewidth',3)
% % p(14)= plot(powellCDsquares(:,1),powellCDsquares(:,2)./1000,'sk','markerfacecolor','k','markersize',8)
% % p(15)= plot(Zhangu10,ZhangCd/1000,'ms','markerfacecolor','m','markersize',8)
% % p(16)= plot(u10plt(:),Andreas_CD(:),'--c','linewidth',4)
% % p(17)=errorbar(HolthuijsenCD(:,1),HolthuijsenCD(:,2),HolthuijsenCD_error,'rs','markerfacecolor','r','markersize',8)
% % 
% % xlabel('\it $U_{10}(m/s)$','FontName','Times New Roman','fontsize',24);
% % ylabel('\it $C_D$','FontName','Times New Roman','fontsize',24);
% % set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',20);
% % set(0,'DefaultTextInterpreter', 'latex');
% % set(gca,'XLim',[0 70]);
% % set(gca,'YLim',[0 0.005]);
% % set(gca,'XTick',[0:10:70]);
% % set(gca,'YTick',[0:0.001:0.005]);
% % 
% % %h1=legend([p(1:6)],'$Z=100$','$Z=150$','$Z=200$','$Z=250$','$Z=300$','$Z=350$');
% % h1=legend([p(1:4)],'$0<R/RMS \le 1$','$1<R/RMS\le2$','$2<R/RMS\le3$','$R/RMS>3$');
% % ah=axes('position',get(gca,'position'),'visible','off');
% % set(h1,'Interpreter','latex','FontSize',16)
% % set(h1,'Box','off');
% % 
% % h2=legend(ah, [p(13:17)],'$Large \& Pond (1981)$','$Powell~et~al.~(2003)$','$Zhang~et~al.~(2008)$','$Andreas~et~al.~(2012)$','$Holthuijsen~et~al.~(2012)$');
% % set(h2,'Interpreter','latex','FontSize',16)
% % set(h2,'Box','off');
% % 
% % eval(['u10_MBL' num2str(Fitting_height) '= u10']); % Copy to new variable with new and different name.
% % clear('u10');  % Delete old variable with old name.
% % eval(['CD_MBL' num2str(Fitting_height) '= CD']);
% % clear('CD');  % Delete old variable with old name.
% % eval(['ustar_MBL' num2str(Fitting_height) '= ustar']);
% % clear('ustar');

end
end
%% Based on above fitting, plot the mean velocity profile for different R/RMW and windbin
% Name: Figures/U_height_RMW_fitting.png
figure
for i=1:1:4
    subplot(2,2,i)
    for j=1:1:8
        p(j)=plot(mean_U_profiles(1:1:end,i,j),zplot(1:1:end),'+')
        hold on
        tmp=mean_U_profiles(:,i,j);
        zfit = zplot(~isnan(tmp));
        U_fit=polyval(Ucoeffs(:,i,j),log(zfit));
        plot(U_fit,zfit,'k-','linewidth',1.5)
%         hold on
%         numtmp = numvecU(:,i,j);
%         stdtmp = std_U_profiles(:,i,j);
%         errorbar(tmp(numtmp>min_samples),zplot(numtmp>min_samples),2*stdtmp(numtmp>min_samples),'*','horizontal')
%         hold on
    end
    
    set(p(1),'linewidth',1,'color','k','linestyle','none');
    set(p(2),'linewidth',1,'color','b','linestyle','none');
    set(p(3),'linewidth',1,'color','r','linestyle','none');
    set(p(4),'linewidth',1,'color','g','linestyle','none');
    set(p(5),'linewidth',1,'color','c','linestyle','none');
    set(p(6),'linewidth',1,'color','k','linestyle','none');
    set(p(7),'linewidth',1,'color','b','linestyle','none');
    set(p(8),'linewidth',1,'color','r','linestyle','none');
    xlabel('\it $wind~speed (m/s)$','FontName','Times New Roman','fontsize',24);
    ylabel('\it $height (m)$','FontName','Times New Roman','fontsize',24);
    set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',20);
    set(0,'DefaultTextInterpreter', 'latex');
    set(gca,'XLim',[0 90]);
    set(gca,'XTick',[0:10:90]);
    set(gca,'YLim',[10 2000]);
    set(gca, 'YScale', 'log')
    box off
end
%% plot the wind profile of different R/RMW in the same wind bin 
% Name: Figures/speed_comparison.png
figure
count=0;
for j=2:1:5
    count=count+1;
    subplot(2,2,count)
    for i=1:1:4
        p(i)=plot(mean_U_profiles(1:1:end,i,j),zplot(1:1:end))
        hold on
    end
    set(p(1),'linewidth',2,'color','k','linestyle','-');
    set(p(2),'linewidth',2,'color','b','linestyle','--');
    set(p(3),'linewidth',2,'color','r','linestyle','-.');
    set(p(4),'linewidth',2,'color','g','linestyle',':');
    xlabel('\it $wind~speed (m/s)$','FontName','Times New Roman','fontsize',24);
    ylabel('\it $height (m)$','FontName','Times New Roman','fontsize',24);
    set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',20);
    set(0,'DefaultTextInterpreter', 'latex');
    %set(gca,'XLim',[0 90]);
    %set(gca,'XTick',[0:10:90]);
    set(gca,'YLim',[10 2000]);
    set(gca, 'YScale', 'log')
    box off
end
    h1=legend(p(1:4),'$0<R/RMW<1$','$1<R/RMW<2$','$2<R/RMW<3$','$R/RMW>3$');
    ah=axes('position',get(gca,'position'),'visible','off');
%     h2=legend(ah,p(5),'$single-phase,Re_b=380$');
%     ah=axes('position',get(gca,'position'),'visible','off');
    set(h1,'Interpreter','latex','FontSize',20)
    set(h1,'Box','off');

%% plot the error of the fitting
% Name: Figures/RMSE_height_RMW_fitting.png
figure
height_fit=[75,100,125,150,175,200,225,250,275,300,325];
for i=1:1:4
    subplot(2,2,i)
    for j=1:1:8
        hold on
        p(j)=plot(height_fit,squeeze(err_fit(i,j,:)),'-o'); 
    end
    set(p(1),'linewidth',1.5,'color','k','linestyle','-');
    set(p(2),'linewidth',1.5,'color','b','linestyle','-');
    set(p(3),'linewidth',1.5,'color','r','linestyle','-');
    set(p(4),'linewidth',1.5,'color','g','linestyle','-');
    set(p(5),'linewidth',1.5,'color','c','linestyle','-');
    set(p(6),'linewidth',1.5,'color','k','linestyle','--');
    set(p(7),'linewidth',1.5,'color','b','linestyle','--');
    set(p(8),'linewidth',1.5,'color','r','linestyle','--');
    xlabel('\it $Fitting~height(m)$','FontName','Times New Roman','fontsize',24);
    ylabel('\it $Root~mean~square~error$','FontName','Times New Roman','fontsize',24);
    set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',20);
    set(0,'DefaultTextInterpreter', 'latex');
    set(gca,'XLim',[50 350]);
    set(gca,'YLim',[0 1.0]);
    h1=legend(p(1:8),'$0-10$','$10-20$','$20-30$','$30-40$','$40-50$','$50-60$','$60-70$','$70-85$');
    ah=axes('position',get(gca,'position'),'visible','off');
%     h2=legend(ah,p(5),'$single-phase,Re_b=380$');
%     ah=axes('position',get(gca,'position'),'visible','off');
    set(h1,'Interpreter','latex','FontSize',20)
    set(h1,'Box','off');
end
%% Based on Dan's data, plot the RMW and corresponding frequency of wind speed.
% Name: Figures/RMW_Rmms.png
load('./RWS/ebtrk_atlc_1988_2017.mat','Year_EBT','Month_EBT','Day_EBT','HourUTC_EBT','StormName_EBT','rmkm_EBT','Vmms_EBT');
V_bin=[0:10:80];
V_center=(V_bin(1:end-1)+V_bin(2:end))/2;
RMW_bin=[0:10:200];
RMW_center=(RMW_bin(1:end-1)+RMW_bin(2:end))/2;
[RMW_V_mesh_RMW,RMW_V_mesh_V]=meshgrid(V_bin,RMW_bin);
RMW_V=[rmkm_EBT(rmkm_EBT>0&Vmms_EBT>0) Vmms_EBT(rmkm_EBT>0&Vmms_EBT>0)];
N = hist3(RMW_V,'Ctrs',{RMW_bin V_bin});
N_sum_RMW=sum(N,1);
N_sum_V=sum(N,2);
figure
subplot(2,2,1)
surf(RMW_V_mesh_V,RMW_V_mesh_RMW,N)
view(2)
colormap(jet)
xlabel('\it $RMW(km)$','FontName','Times New Roman','fontsize',24);
ylabel('\it $Max~velocity(m/s)$','FontName','Times New Roman','fontsize',24);
set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',20);
set(0,'DefaultTextInterpreter', 'latex');
subplot(2,2,3)
plot(V_center,N_sum_RMW(2:end),'k*','markersize',8)
xlabel('\it $Max~velocity(m/s)$','FontName','Times New Roman','fontsize',24);
ylabel('\it $Frequency$','FontName','Times New Roman','fontsize',24);
set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',20);
set(0,'DefaultTextInterpreter', 'latex');
subplot(2,2,4)
plot(RMW_center,N_sum_V(2:end),'k*','markersize',8)
xlabel('\it $RMW(km)$','FontName','Times New Roman','fontsize',24);
ylabel('\it $Frequency$','FontName','Times New Roman','fontsize',24);
set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',20);
set(0,'DefaultTextInterpreter', 'latex');
% view(2)
%% Based on Dan's data, choose three hurricanes, plot their RMW(t)
% Name: Figures/EPSILON_HARVEY_MARIA_RMW_VS_time.png
load('./RWS/ebtrk_atlc_1988_2017.mat','Year_EBT','Month_EBT','Day_EBT','HourUTC_EBT','StormName_EBT','rmkm_EBT','Vmms_EBT');
A=rmkm_EBT(7532:7574); %EPSILON2005
B=rmkm_EBT(12359:12428);%HARVEY2017
C=rmkm_EBT(12662:12726);%MARIA2017
figure
plot(A,'k-+','linewidth',1.5)
hold on
plot(B,'b-*','linewidth',1.5)
hold on
plot(C,'r-o','linewidth',1.5)
xlabel('\it $Every~6~hours$','FontName','Times New Roman','fontsize',24);
ylabel('\it $RMW(km)$','FontName','Times New Roman','fontsize',24);
set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',20);
set(0,'DefaultTextInterpreter', 'latex');
h1=legend('$EPSILON2005$','$HARVEY2017$','$MARIA2017$');
ah=axes('position',get(gca,'position'),'visible','off');
set(h1,'Interpreter','latex','FontSize',14)
set(h1,'Box','off');
%% Based on output data, CD_ustar_MBL100-500.mat, plot u* based on different fitting height.
% Name: Figures/ustar_RMW_MBL100-500.png,(a) R/RMW<1, (b) 1<R/RMW<2, (c)2<R/RMW<3, (d) R/RMW>3
load ./CD_ustar_MBL100-500.mat
powellCDsquares = load('./Published_data/powellCDsquares.dat');
Powellustar = load('./Published_data/Powell_ustar.dat');
figure
for i=1:1:4
    subplot(2,2,i)
    p(1)=plot(u10_MBL100(i,:),ustar_MBL100(i,:),'k-*','linewidth',1.5);
    hold on
    p(2)=plot(u10_MBL200(i,:),ustar_MBL200(i,:),'b-*','linewidth',1.5);
    hold on
    p(3)=plot(u10_MBL300(i,:),ustar_MBL300(i,:),'r-*','linewidth',1.5);
    hold on
    p(4)=plot(u10_MBL400(i,:),ustar_MBL400(i,:),'g-*','linewidth',1.5);
    hold on
    p(5)=plot(u10_MBL500(i,:),ustar_MBL500(i,:),'m-*','linewidth',1.5);
    hold on
    p(6)= plot(Powellustar(:,1),Powellustar(:,2),'sk','markerfacecolor','k','markersize',8)
xlabel('\it $U_{10}(m/s)$','FontName','Times New Roman','fontsize',24);
ylabel('\it $Friction~velocity(m~s^{-1})$','FontName','Times New Roman','fontsize',24);
set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',20);
set(0,'DefaultTextInterpreter', 'latex');
set(gca,'XLim',[0 70]);
set(gca,'YLim',[0 3]);
set(gca,'XTick',[0:10:70]);
set(gca,'YTick',[0:1.0:3.0]);
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickLength',[0.015, 0.015])
end
h1=legend(p(1:6),'$MBL:100M,Fitting:100m$','$MBL:200M,Fitting:200m$','$MBL:300M,Fitting:300m$'...
                ,'$MBL:400M,Fitting:400m$','$MBL:500M,Fitting:500m$','$Powell~et~al,2003$');
ah=axes('position',get(gca,'position'),'visible','off');
set(h1,'Interpreter','latex','FontSize',14)
set(h1,'Box','off');
%% Based on output data, CD_MBL100-500.mat, plot Cd based on different fitting height.
% Name: Figures/Cd_RMW_MBL100-500.png, (a) R/RMW<1, (b) 1<R/RMW<2, (c)2<R/RMW<3, (d) R/RMW>3
load ./CD_MBL100-500.mat
powellCDsquares = load('./Published_data/powellCDsquares.dat');
Powellustar = load('./Published_data/Powell_ustar.dat');
figure
for i=1:1:4
    subplot(2,2,i)
    p(1)=plot(u10_MBL100(i,:),CD_MBL100(i,:),'k-*','linewidth',1.5);
    hold on
    p(2)=plot(u10_MBL200(i,:),CD_MBL200(i,:),'b-*','linewidth',1.5);
    hold on
    p(3)=plot(u10_MBL300(i,:),CD_MBL300(i,:),'r-*','linewidth',1.5);
    hold on
    p(4)=plot(u10_MBL400(i,:),CD_MBL400(i,:),'g-*','linewidth',1.5);
    hold on
    p(5)=plot(u10_MBL500(i,:),CD_MBL500(i,:),'m-*','linewidth',1.5);
    hold on
    p(6)= plot(powellCDsquares(:,1),powellCDsquares(:,2)./1000,'sk','markerfacecolor','k','markersize',8)
xlabel('\it $U_{10}(m/s)$','FontName','Times New Roman','fontsize',24);
ylabel('\it $C_D$','FontName','Times New Roman','fontsize',24);
set(gca,'FontName','Times New Roman','linewidth',1.5,'fontsize',20);
set(0,'DefaultTextInterpreter', 'latex');
set(gca,'XLim',[0 70]);
set(gca,'YLim',[0 0.003]);
set(gca,'XTick',[0:10:70]);
set(gca,'YTick',[0:0.001:0.003]);
end
h1=legend(p(1:6),'$MBL:100M,Fitting:100m$','$MBL:200M,Fitting:200m$','$MBL:300M,Fitting:300m$'...
                ,'$MBL:400M,Fitting:400m$','$MBL:500M,Fitting:500m$','$Powell~et~al,2003$');
ah=axes('position',get(gca,'position'),'visible','off');
set(h1,'Interpreter','latex','FontSize',14)
set(h1,'Box','off');


%%
figure(5)
for i=1:num_rad_bins
    tmp=mean_U_profiles(:,i);
    stdtmp = std_U_profiles(:,i);
    zfit = zplot(~isnan(tmp));
    subplot(1,num_rad_bins,i)
    hold all
    numtmp = numvecU(:,i);
    
    errorbar(zplot(numtmp>min_samples),tmp(numtmp>min_samples),2*stdtmp(numtmp>min_samples),'*')
    plot(zfit,polyval(Ucoeffs(:,i),log(zfit)),'-')
    set(gca,'xscale','log')
    ylabel('$<U> (m/s)$')
    xlabel('$z (m)$')
    
%     fnum1 = fopen(['./dat_files/umean' num2str(i) '.dat'],'w');
%     fnum2 = fopen(['./dat_files/umean' num2str(i) '_fit.dat'],'w');
%     for k=1:length(zplot(numtmp>min_samples))
%         if (numtmp(i)>min_samples)
%             fprintf(fnum1,'%e\t%e\t%e\n',zplot(k),tmp(k),2*stdtmp(k));
%         end
%     end
%     ufit_out = polyval(Ucoeffs(:,i),log(zfit));
%     for k=1:length(zfit)
%         fprintf(fnum2,'%e\t%e\n',zfit(k),ufit_out(k));
%     end
%     fclose(fnum1);
%     fclose(fnum2);    
end
%% III save pic
% save fig, jpeg, png, pdf in folder ./Figures
save=1;
%filename='./../FIGURES';%filename
filename='./Figures';%filename
strFileName = 'ustar_RMW_MBL100-500';%figure name
paper_size='A3';
paper_direction='landscape'; %'portrait'
    % Fragment: save figure with fig png jpg
    Frag_figure_save