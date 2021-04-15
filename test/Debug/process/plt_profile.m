clc;clear;

fdir1='../OP2/';
GISfile='OP2.gif';

dep=load('../Input/dep2.txt');

s2g=load([fdir1,'sta_grn_0002']);
s1=load([fdir1,'sta_0001']);

m=66;

ll=[-0.07:0.01:0.07];
len=length(ll);

xd=[0.25:0.25:m];
xdForFill=[xd,fliplr(xd)];
botd=zeros(size(xd))-0.214;
ee=-dep(12,:);
ydForFill=[botd,fliplr(ee)];


x=[1:m];
xForFill=[x,fliplr(x)];
bot=zeros(size(x))-0.214;

stepall=71;
for i=1:stepall
% i=355;
    j=i+99;
%     eta=load([fdir1,'eta_0',num2str(j)]);
    etag=load([fdir1,'etagrn_0',num2str(j)]);
    
%     ee(1:m/2)=eta(2,1:m/2);
    ee=etag(2,:);
    
    tt=zeros(len,1)+10*j;
    
    subplot(211)
    plot(s1(:,1),s1(:,2),'LineWidth',2);
    hold on
    plot(s2g(:,1),s2g(:,2),'r--','LineWidth',2);
    hold on
    plot(tt,ll,'k','LineWidth',2)
    
    xlim([950 1810])
    ylim([-0.07 0.07])
    ylabel('Waterlevel(m)')
    title(['t = ',num2str(j*10),'s'])
    
    subplot(212)
    fill(xdForFill,ydForFill,'k','FaceAlpha',0.5,'EdgeAlpha',0.5,'EdgeColor','k');
    hold on
    yForFill=[ee,bot];
    fill(xForFill,yForFill,'c','FaceAlpha',0.5,'EdgeAlpha',0.5,'EdgeColor','c');
    xlim([1 m])
    ylim([-0.2 0.1])
    
    ylabel('Waterlevel\eta(m)');
    
    set(gcf,'color','w');
    picname=strcat(num2str(j),'.fig');
    picname=char(picname);
    saveas(gcf,picname);
    
    frame=getframe(gcf);  
    im=frame2im(frame);
    [I,map]=rgb2ind(im,20);          
    if i==1
        imwrite(I,map,GISfile,'gif', 'Loopcount',inf,'DelayTime',0.5);
    elseif i==stepall
        imwrite(I,map,GISfile,'gif','WriteMode','append','DelayTime',0.5);
    else
        imwrite(I,map,GISfile,'gif','WriteMode','append','DelayTime',0.5);

    end;  
    close all
      
end











        
        
        