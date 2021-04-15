clear all
dx=0.02;

x=[-2.64:dx:2.64-dx]; % 264
% y=[-0.02:0.02:0.02];

m=length(x);
% n=length(y);
n=0.48/dx;


% s=0.214+min(0,max(-0.33,0.5*abs(x-1.4)-0.35));
s=0.214+min(0,max(-0.250,0.5*abs(x-1.4)-0.35));

dep=zeros(n,m);
for j=1:n
dep(j,:)=s(:);
end

figure(1)
plot(x,-dep(2,:))
figure(2)
pcolor(dep)
shading interp

save -ASCII dep2.txt dep