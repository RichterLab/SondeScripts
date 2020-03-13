% This code loads in trajectory data from CM1

tic;

% Filename containing trajectory data
filename=['./cm1out_pdata.nc'];

Rd=287.04;
cp=1005.7;

time_traj=ncread(filename,'time');
xh=ncread(filename,'xh');
yh=ncread(filename,'yh');
zh=ncread(filename,'zh');

x=ncread(filename,'x');
y=ncread(filename,'y');
z=ncread(filename,'z');
u=ncread(filename,'u');
v=ncread(filename,'v');
w=ncread(filename,'w');
theta=ncread(filename,'th');
pres=ncread(filename,'prs');
b=ncread(filename,'b');
dpdz=ncread(filename,'dpdz');
qv=ncread(filename,'qv');
qc=ncread(filename,'qc');
%qr=ncread(filename,'qr');
%qi=ncread(filename,'qi');
%qs=ncread(filename,'qs');
%qg=ncread(filename,'qg');
kh=ncread(filename,'kh');

dbz=ncread(filename,'dbz');

rho=ncread(filename,'rho');

zeta=ncread(filename,'zvort');
xvort=ncread(filename,'xvort');
yvort=ncread(filename,'yvort');

zvort_stretching=ncread(filename,'zv1');
zvort_tilting=ncread(filename,'zv2');
zvort_nonbous=ncread(filename,'zv3');
zvort_turbdif=ncread(filename,'zv4');
zvort_impdif=ncread(filename,'zv5');

wspd=sqrt(u.^2+v.^2);

numtraj=size(u,1);
numtimes=size(u,2);
for j=1:numtraj
if min(z(j,:))<=1
ind_end(j)=find(wspd(j,:)<=1e-4,1);
x(j,ind_end(j):end)=NaN;
y(j,ind_end(j):end)=NaN;
z(j,ind_end(j):end)=NaN;
u(j,ind_end(j):end)=NaN;
v(j,ind_end(j):end)=NaN;
wspd(j,ind_end(j):end)=NaN;
w(j,ind_end(j):end)=NaN;
theta(j,ind_end(j):end)=NaN;
pres(j,ind_end(j):end)=NaN;
b(j,ind_end(j):end)=NaN;
dpdz(j,ind_end(j):end)=NaN;
qv(j,ind_end(j):end)=NaN;
qc(j,ind_end(j):end)=NaN;
%qr(j,ind_end(j):end)=NaN;
kh(j,ind_end(j):end)=NaN;
dbz(j,ind_end(j):end)=NaN;
rho(j,ind_end(j):end)=NaN;
zeta(j,ind_end(j):end)=NaN;
xvort(j,ind_end(j):end)=NaN;
yvort(j,ind_end(j):end)=NaN;
zvort_stretching(j,ind_end(j):end)=NaN;
zvort_tilting(j,ind_end(j):end)=NaN;
zvort_nonbous(j,ind_end(j):end)=NaN;
zvort_turbdif(j,ind_end(j):end)=NaN;
zvort_impdif(j,ind_end(j):end)=NaN;

time_to_sfc(j)=time_traj(ind_end(j))-time_traj(1);

else 
ind_end(j)=201;
time_to_sfc(j)=3+time_traj(end)-time_traj(1);

%ind_z10=find(abs(z(j,:)-10)==min(abs(z(j,:)-10)),1);
%z10(j)=z(j,ind_z10);
%wspd10(j)=wspd(j,ind_z10);

end
end


temp=theta.*((pres./100000).^(Rd/cp));


%wspd=sqrt(u.^2+v.^2);


%time_center=ncread(filename_center,'time');
time_center=[0:60:21240];
%icenter=squeeze(ncread(filename_center,'icenter'));
%jcenter=squeeze(ncread(filename_center,'jcenter'));
%xcenter=squeeze(ncread(filename_center,'xcenter'));
%ycenter=squeeze(ncread(filename_center,'ycenter'));

%vbarmax=squeeze(ncread(filename_center,'vmax'));
%rbarmax=squeeze(ncread(filename_center,'rmax'));
%zbarmax=squeeze(ncread(filename_center,'zmax'));

%icenter(end)=[];
%jcenter(end)=[];
%xcenter(end)=[];
%ycenter(end)=[];
%vbarmax(end)=[];
%rbarmax(end)=[];
%zbarmax(end)=[];

xcenter=zeros(size(time_center));
ycenter=zeros(size(time_center));

xcenter_interp=interp1(time_center,xcenter,time_traj);
ycenter_interp=interp1(time_center,ycenter,time_traj);
for j=1:length(time_traj)
 xdist(:,j)=x(:,j)-xcenter(j);
 ydist(:,j)=y(:,j)-ycenter(j);
 alpha(:,j)=atan2(ydist(:,j),xdist(:,j));
 vt(:,j)=cos(alpha(:,j)).*v(:,j)-sin(alpha(:,j)).*u(:,j);
 vr(:,j)=cos(alpha(:,j)).*u(:,j)+sin(alpha(:,j)).*v(:,j);
end

dist=sqrt(xdist.^2+ydist.^2);

%dist=sqrt(x.^2+y.^2);

[r5to15]=find(dist(:,1)<=15000 & dist(:,1)>=5000);
[r0to5]=find(dist(:,1)<=5000);
[r5to10]=find(dist(:,1)<=10000 & dist(:,1)>=5000);
[r10to15]=find(dist(:,1)<=15000 & dist(:,1)>=10000);
[r15to20]=find(dist(:,1)<=20000 & dist(:,1)>=15000);
[r20to25]=find(dist(:,1)<=25000 & dist(:,1)>=20000);
[r25to30]=find(dist(:,1)<=30000 & dist(:,1)>=25000);

num_w10_r0to5=size(find(max(w(r0to5,:),[],2)>=10));
num_w10_r5to10=size(find(max(w(r5to10,:),[],2)>=10));
num_w10_r10to15=size(find(max(w(r10to15,:),[],2)>=10));
num_w10_r15to20=size(find(max(w(r15to20,:),[],2)>=10));
num_w10_r20to25=size(find(max(w(r20to25,:),[],2)>=10));
num_w10_r25to30=size(find(max(w(r25to30,:),[],2)>=10));

pct_w10_r0to5=num_w10_r0to5(1)/length(r0to5);
pct_w10_r5to10=num_w10_r5to10(1)/length(r5to10);
pct_w10_r10to15=num_w10_r10to15(1)/length(r10to15);
pct_w10_r15to20=num_w10_r15to20(1)/length(r15to20);
pct_w10_r20to25=num_w10_r20to25(1)/length(r20to25);
pct_w10_r25to30=num_w10_r25to30(1)/length(r25to30);

%num_w10to20_r0to5=size(find(max(w(r0to5,:),[],2)>=10));
%num_w10to20_r5to10=size(find(max(w(r5to10,:),[],2)>=10));
%num_w10to20_r10to15=size(find(max(w(r10to15,:),[],2)>=10 & max(w(r10to15,:),[],2)<20));
%num_w10_r15to20=size(find(max(w(r15to20,:),[],2)>=10));
%num_w10_r20to25=size(find(max(w(r20to25,:),[],2)>=10));
%num_w10_r25to30=size(find(max(w(r25to30,:),[],2)>=10));

%pct_w10_r0to5=num_w10_r0to5(1)/length(r0to5);
%pct_w10_r5to10=num_w10_r5to10(1)/length(r5to10);
%pct_w10_r10to15=num_w10_r10to15(1)/length(r10to15);
%pct_w10_r15to20=num_w10_r15to20(1)/length(r15to20);
%pct_w10_r20to25=num_w10_r20to25(1)/length(r20to25);
%pct_w10_r25to30=num_w10_r25to30(1)/length(r25to30);

num_wspd100_r0to5=size(find(max(wspd(r0to5,:),[],2)>=100));
num_wspd100_r5to10=size(find(max(wspd(r5to10,:),[],2)>=100));
num_wspd100_r10to15=size(find(max(wspd(r10to15,:),[],2)>=100));
num_wspd100_r15to20=size(find(max(wspd(r15to20,:),[],2)>=100));
num_wspd100_r20to25=size(find(max(wspd(r20to25,:),[],2)>=100));
num_wspd100_r25to30=size(find(max(wspd(r25to30,:),[],2)>=100));

pct_wspd100_r0to5=num_wspd100_r0to5(1)/length(r0to5);
pct_wspd100_r5to10=num_wspd100_r5to10(1)/length(r5to10);
pct_wspd100_r10to15=num_wspd100_r10to15(1)/length(r10to15);
pct_wspd100_r15to20=num_wspd100_r15to20(1)/length(r15to20);
pct_wspd100_r20to25=num_wspd100_r20to25(1)/length(r20to25);
pct_wspd100_r25to30=num_wspd100_r25to30(1)/length(r25to30);


%z_interp=[0:10:3000];
%for j=1:numtraj
% w_interp=interp1(z(j,1:ind_end(j)-1),w(j,1:ind_end(j)-1),z_interp);
%end

save traj_test x y z

toc;
