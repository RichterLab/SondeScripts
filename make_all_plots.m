
clc; close all

%Must run Bin_SondeData to generate the following files:

load('constants.mat');
load('all_sst.mat');
load('all_profile_data.mat');

min_samples = 2;   %Minimum number of samples to be included
min_fit_samples = 5;

zplot = min_height+0.5*height_interval:height_interval:max_height-0.5*height_interval;

%For fitting the log region in these mean profiles:
Ucoeffs = zeros(2,num_wind_bins,num_rad_bins);
thetacoeffs = zeros(2,num_wind_bins,num_theta_bins,num_rad_bins);
qcoeffs = zeros(2,num_wind_bins,num_theta_bins,num_rad_bins);
kcoeffs = zeros(2,num_wind_bins,num_theta_bins,num_rad_bins);

ustar = zeros(num_wind_bins,num_rad_bins);
u10 = zeros(num_wind_bins,num_rad_bins);
CD = zeros(num_wind_bins,num_rad_bins);
thetastar = zeros(num_wind_bins,num_theta_bins,num_rad_bins);
CH = zeros(num_wind_bins,num_theta_bins,num_rad_bins);
qstar = zeros(num_wind_bins,num_theta_bins,num_rad_bins);
CE = zeros(num_wind_bins,num_theta_bins,num_rad_bins);
kstar = zeros(num_wind_bins,num_theta_bins,num_rad_bins);
CK = zeros(num_wind_bins,num_theta_bins,num_rad_bins);

mean_U_profiles = zeros(num_z,num_wind_bins,num_rad_bins);
std_U_profiles = zeros(num_z,num_wind_bins,num_rad_bins);
numvecU = zeros(num_z,num_wind_bins,num_rad_bins);
mean_theta_profiles = zeros(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
std_theta_profiles = zeros(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
numvecT = zeros(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
mean_q_profiles = zeros(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
std_q_profiles = zeros(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
mean_k_profiles = zeros(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
std_k_profiles = zeros(num_z,num_wind_bins,num_theta_bins,num_rad_bins);
mean_p_profiles = zeros(num_z,num_wind_bins,num_theta_bins,num_rad_bins);

% Setting up data to plot CD and CK lines
U10 = linspace(0,80,100);

%These are the actual surface flux coefficients in cm1:
dcd1 = 1.0e-3;
dcd2 = 2.4e-3;
dwsp1 = 5.0;
dwsp2 = 25.0;

for i=1:length(U10)
    cd(i) = dcd1 + (U10(i) - dwsp1)*(dcd2-dcd1)/(dwsp2-dwsp1);
    cd(i) = min(cd(i),dcd2);
    cd(i) = max(cd(i),dcd1);
    ck(i) = 1.2e-3;
end


for i=1:num_wind_bins
    for r=1:num_rad_bins
        
        for k=1:num_z
            mean_U_profiles(k,i,r) = mean(all_U_profiles{k,i,r}(:));
            std_U_profiles(k,i,r) = std(all_U_profiles{k,i,r}(:));
            numvecU(k,i,r) = length(all_U_profiles{k,i,r}(:));
        end
        
        for j = 1:num_theta_bins
            for k=1:num_z
                mean_theta_profiles(k,i,j,r) = mean(all_theta_profiles{k,i,j,r}(:));
                std_theta_profiles(k,i,j,r) = std(all_theta_profiles{k,i,j,r}(:));
                numvecT(k,i,j,r) = length(all_theta_profiles{k,i,j,r}(:));
                
                mean_q_profiles(k,i,j,r) = mean(all_q_profiles{k,i,j,r}(:));
                std_q_profiles(k,i,j,r) = std(all_q_profiles{k,i,j,r}(:));
                
                mean_k_profiles(k,i,j,r) = mean(all_k_profiles{k,i,j,r}(:));
                std_k_profiles(k,i,j,r) = std(all_k_profiles{k,i,j,r}(:));
                
                mean_p_profiles(k,i,j,r) = mean(all_p_profiles{k,i,j,r}(:));
            end
        end
    end
end


for i=1:num_wind_bins
    for r=1:num_rad_bins
        
        tmp = mean_U_profiles(:,i,r);
        numtmp = numvecU(:,i,r);
        zfit = zplot(~isnan(tmp)&numtmp>min_samples);
        ufit = tmp(~isnan(tmp)&numtmp>min_samples);
        if (length(ufit) > min_fit_samples)
            Ucoeffs(:,i,r) = polyfit(log(zfit),ufit',1);
            test = fit(log(zfit)',ufit,'poly1');
            U_ci = confint(test,0.95);
        else
            Ucoeffs(:,i,r) = ones(2,1)*NaN;
            U_ci = ones(2,1)*NaN;
        end
        
        
        
        %Compute the mean velocity at the 10m using fit:
        u10(i,r) = Ucoeffs(1,i,r)*log(10) + Ucoeffs(2,i,r);
        
        ustar(i,r) = Ucoeffs(1,i,r)*0.4;
        CD(i,r) = ustar(i,r)^2/u10(i,r)^2;
        
        %compute errors:
        delta_u10(i,r) = 2*std_U_profiles(1,i,r);
        delta_ustar(i,r) = 0.5*(U_ci(2,1) - U_ci(1,1));
        delta_CD(i,r) = CD(i,r)*sqrt(2*(delta_ustar(i,r)/abs(ustar(i,r)))^2 + 2*(delta_u10(i,r)/abs(u10(i,r)))^2);
        
        for j=1:num_theta_bins
            
            %Fit the theta profile:
            numtmp = numvecT(:,i,j,r);
            thetatmp = mean_theta_profiles(:,i,j,r);
            zthetafit = zplot(~isnan(thetatmp)&numtmp>min_samples);
            thetafit = thetatmp(~isnan(thetatmp)&numtmp>min_samples);
            if (length(thetafit) > min_fit_samples)
                thetacoeffs(:,i,j,r) = polyfit(log(zthetafit),thetafit',1);
                test = fit(log(zthetafit)',thetafit,'poly1');
                theta_ci = confint(test,0.95);
            else
                thetacoeffs(:,i,j,r) = ones(2,1)*NaN;
                theta_ci = zeros(2,2);
            end

            theta10 = thetacoeffs(1,i,j,r)*log(10) + thetacoeffs(2,i,j);
            
            theta_s = mean(ssttosave{i,j,r});
            
            thetastar(i,j,r) = thetacoeffs(1,i,j,r)*0.4;
            CH(i,j,r) = -ustar(i,r)*thetastar(i,j,r)/u10(i,r)/(theta_s-theta10);
            
            %compute errors:
            delta_theta10 = 2*std_theta_profiles(1,i,j,r);
            delta_thetastar(i,j,r) = 0.5*(theta_ci(2,1) - theta_ci(1,1));
            delta_sst = 2*std(ssttosave{i,j,r});
            Deltatheta = (theta_s - theta10);
            delta_Deltatheta = sqrt(delta_sst^2 + delta_theta10^2);
            delta_CH(i,j,r) = CH(i,j,r)*sqrt((delta_ustar(i,r)/abs(ustar(i,r)))^2 + ...
                (delta_thetastar(i,j,r)/abs(thetastar(i,j,r)))^2 + ...
                (delta_u10(i,r)/abs(u10(i,r)))^2 + ...
                (delta_Deltatheta/abs(Deltatheta))^2);
            
            %Fit the moisture profile:
            numtmp = numvecT(:,i,j,r);
            qtmp = mean_q_profiles(:,i,j,r);
            zqfit = zplot(~isnan(qtmp)&numtmp>min_samples);
            qfit = qtmp(~isnan(qtmp)&numtmp>min_samples);
            if (length(qfit) > min_fit_samples)
                qcoeffs(:,i,j,r) = polyfit(log(zqfit),qfit',1);
                test = fit(log(zqfit)',qfit,'poly1');
                q_ci = confint(test,0.95);
            else
                qcoeffs(:,i,j,r) = ones(2,1)*NaN;
                q_ci = zeros(2,2);
            end
            
            q10 = qcoeffs(1,i,j,r)*log(10) + qcoeffs(2,i,j,r);
            
            %Compute the surface moisture (assuming saturation) based on the
            %surface temperature theta_s:
            estar_s = 6.11*exp(17.3*(theta_s-273)./((theta_s-273)+237.3));
            
            psurf = mean_p_profiles(1,i,j);
            eoverp = estar_s/psurf;
            q_s = 0.622*eoverp/(1-0.378*eoverp)*100;
            %q_s = qcoeffs(2,i,j);
            qstar(i,j,r) = qcoeffs(1,i,j,r)*0.4;
            CE(i,j,r) = -ustar(i,r)*qstar(i,j,r)/u10(i,r)/(q_s-q10);
            
            %Need to convert SST to qs and ks to get a distribution:
            estar_s_vec = 6.11*exp(17.3*(ssttosave{i,j,r}-273)./((ssttosave{i,j,r}-273)+237.3));
            eoverp_vec = estar_s_vec/psurf;
            qs_vec = 0.622*eoverp_vec./(1-0.378*eoverp);
            qs_std = std(qs_vec);
            
            %compute errors:
            delta_q10 = 2*std_q_profiles(1,i,j,r);
            delta_qstar(i,j,r) = 0.5*(q_ci(2,1) - q_ci(1,1));
            delta_qs = 2*qs_std;
            Deltaq = (q_s - q10);
            delta_Deltaq = sqrt(delta_qs^2 + delta_q10^2);
            delta_CE(i,j,r) = CE(i,j,r)*sqrt((delta_ustar(i,r)/abs(ustar(i,r)))^2 + ...
                (delta_qstar(i,j,r)/abs(qstar(i,j,r)))^2 + ...
                (delta_u10(i,r)/abs(u10(i,r)))^2 + ...
                (delta_Deltaq/abs(Deltaq))^2);
            
            %Fit the enthalpy profile:
            numtmp = numvecT(:,i,j,r);
            ktmp = mean_k_profiles(:,i,j,r);
            zkfit = zplot(~isnan(ktmp)&numtmp>min_samples);
            kfit = ktmp(~isnan(ktmp)&numtmp>min_samples);

            if (length(kfit) > min_fit_samples)
                kcoeffs(:,i,j,r) = polyfit(log(zkfit),kfit',1);
                test = fit(log(zkfit)',kfit,'poly1');
                k_ci = confint(test,0.95);
            else
                kcoeffs(:,i,j,r) = ones(2,1)*NaN;
                k_ci = zeros(2,2);
            end
            
            k10 = kcoeffs(1,i,j)*log(10) + kcoeffs(2,i,j);
            %Compute the surface enthalpy based on surface q, theta:
            k_s = (cpa*(1-q_s) + cpl*q_s)*theta_s + Lv*q_s;
            ks_vec = (cpa*(1-qs_vec) + cpl*qs_vec).*ssttosave{i,j,r} + Lv*qs_vec;
            ks_std = std(ks_vec);
            kstar(i,j,r) = kcoeffs(1,i,j,r)*0.4;
            CK(i,j,r) = -ustar(i,r)*kstar(i,j,r)/u10(i,r)/(k_s-k10);
            
            %compute errors:
            delta_k10 = 2*std_k_profiles(1,i,j,r);
            delta_kstar(i,j,r) = 0.5*(k_ci(2,1) - k_ci(1,1));
            delta_ks = 2*ks_std;
            Deltak(i,j,r) = (k_s - k10);
            delta_Deltak(i,j,r) = sqrt(delta_ks^2 + delta_k10^2);
            ksmean(i,j,r) = mean(ks_vec);
            ksstd(i,j,r) = std(ks_vec);
            k10mean(i,j,r) = mean_k_profiles(1,i,j,r);
            k10std(i,j,r) = std_k_profiles(1,i,j,r);
            
            delta_CK(i,j,r) = CK(i,j,r)*sqrt((delta_ustar(i)/abs(ustar(i)))^2 + ...
                (delta_kstar(i,j,r)/abs(kstar(i,j,r)))^2 + ...
                (delta_u10(i,r)/abs(u10(i,r)))^2 + ...
                (delta_Deltak(i,j,r)/abs(Deltak(i,j,r))^2));
        end
    end
end


figure(1)
h = 1;
i=1;
for r=1:num_rad_bins
        
        tmp=mean_U_profiles(:,i,r);
        stdtmp = std_U_profiles(:,i,r);
        zfit = zplot(~isnan(tmp));
        subplot(1,num_rad_bins,h)
        hold all
        numtmp = numvecU(:,i,r);
        
        plot(tmp(numtmp>min_samples),zplot(numtmp>min_samples),'*')
        plot(polyval(Ucoeffs(:,i,r),log(zfit)),zfit,'-')
        set(gca,'yscale','log')
        xlabel('<U> (m/s)')
        ylabel('z (m)')
        h=h+1;
end



