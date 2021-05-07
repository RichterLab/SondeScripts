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
        meanWS500(h) = nanmean(all_U_profiles(ht500,h));
    else
        meanWS500(h) = NaN;
    end
    numZ(h) = length(ht500);
    clear ht500
end

radRMW = radius_km./RMW;
madhigh = find(meanWS500 >= 70 & minZ <= 150 & radRMW <= 10)


figure(6)
set(gcf,'position',[100 100 800 400])
plotYear = 1996.5:1:2018.5;
histogram(hurrYear(madhigh),plotYear,'facecolor',[0.4 0.9 0.3]);
set(gca,'fontsize',14)
hold on
xarrow = [0.44,0.44];
yarrow = [0.52,0.42];
text(2005,26,{'GPS receiver';'updated'},'rotation',90,'fontsize',14);
annotation('arrow',xarrow,yarrow,'linewidth',2)
xarrow = [0.598,0.598];
yarrow = [0.24,0.14];
text(2010,10,{'wind sampling';'increased from';'2 Hz to 4 Hz'},'rotation',90,'fontsize',14);
annotation('arrow',xarrow,yarrow,'linewidth',2)
xarrow2 = [0.533,0.533];
yarrow2 = [0.22,0.12];
text(2008,8,{'upgrade to';'AVAPS II'},'rotation',90,'fontsize',14);
annotation('arrow',xarrow2,yarrow2,'linewidth',2)
xlabel('Year')
ylabel('Number of sondes with {\it{U}}_{mean} > 70 m s^{-1}')


