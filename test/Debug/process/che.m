fdir1='../OP1/';
fdir2='../OP2/';
sta1=load([fdir2,'sta_grn_0001']);
sta2=load([fdir2,'sta_grn_0002']);

figure
plot(sta1(:,1),sta1(:,2))
hold on
plot(sta2(:,1),sta2(:,2))


% ee=load([fdir2,'etagrn_0360']);
% figure
% pcolor(ee)
% shading flat
