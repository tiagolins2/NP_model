%function f=xDLVO_1(I,x,fthe)
%function f=xDLVO_1(I)
function f=pb_het_DLVO(I,T,qsi, phip, Ii, rho_lig, Mw,Anp,Ass,Rnp,Rss)
x=[1.012706775	0.99708169	1.010199837	1.021996645	1	1	1	1	1	1];

% function f=xDLVO_1(I,a, T,qsi, phip, Ii, rho_lig)
% x=1.2752;
% y=1.157;
%x=[1 1.00844263680481 1.05516790932389 0.939787014129714 1 1 1 1 1 1 1];
%x=[1 1.00844263680481 1.05516790932389 0.939787014129714 1 1 1 1 1 1 1];
A=abs((sqrt(Anp)-sqrt(3.7e-20))*(sqrt(Ass)-sqrt(3.7e-20)));
Kb=1.38e-23;
%a=6.5e-9;
%a=100/2*1e-9;%*x(10);
a=2*Rnp*Rss/(Rnp+Rss);
zp=-3.93e-2;
epsilon0=8.854187e-12;
epsilon=78.4098836*x(1);
%A=3.5e-21*x(2);
%A=1E-20*x(1);
Na=6.0221476E+23;
%T=298.15;
e=1.6e-19;
%I=450;
Z=1;
%qsi=0.45*x(5);
%phip=0.0*x(6);
lambda=1e-7;
 %Ii=1.0e-9; s=5e-8*x(4);
 %Mw=148.22;
% rho_lig=1.087;
%Ii=1.2e-9*x(7); 
s=5e-8; 
%zp=(-0.00000000000037588*I^4+0.0000000007404*I^3-0.0000005107*I^2+0.00014995*I-0.034438)*2.171318648*x(2);
zp=0.001*(7.1356*log(I)*x(7)-62.249*x(8))*1.15;
 if I<20
    zp=-0.0375*x(9);
 end
if Ii==0
   Ii=0.001; phip=0; 
end
h=linspace(0.001,20,10000)*1e-9;
bh=zeros(1,numel(h));
Vtot=zeros(1,numel(h));
Vedl=zeros(1,numel(h));
Vvdw=zeros(1,numel(h));
Vo=zeros(1,numel(h));
Ve=zeros(1,numel(h));
Va=zeros(1,numel(h));
Vr=zeros(1,numel(h));
Fst=zeros(1,numel(h));
Vst=zeros(1,numel(h));
dVtot=zeros(1,numel(h));
int_Vtot=zeros(1,numel(h));
int_Vvdw=zeros(1,numel(h));
dh=1;
Cinv=sqrt((2*I*Na*e^2)/(epsilon*epsilon0*Kb*T));
n_inf=1000*Na*I*0.001;
phi=Z*e*zp/(Kb*T); kappa=2.32e9*(I*0.001*Z^2)^0.5; b=5.32;
kappa=Cinv;
zp=zp*(1+1/(kappa*a));
while dh<=numel(h)
   Vvdw(dh)=A*Rnp*Rss/(6*h(dh)*(Rnp+Rss));
   Va(dh)=A*a/(12*h(dh))*(1-b*h(dh)/lambda*log(1+lambda/(b*h(dh))));
   Vr(dh)=2*pi*a*n_inf*Kb*T*phi^2/kappa^2*(log((1+exp(-kappa*h(dh)))/(1-exp(-kappa*h(dh))))+log(1-exp(-2*kappa*h(dh))));
   gamma=tanh(Z*e*zp/(4*Kb*T));
   
   Vedl(dh)=64*pi*epsilon0*epsilon*Kb^2*T^2*(Rnp*Rss/(Rnp+Rss))/(Z^2*e^2)*gamma^2*exp(-Cinv*h(dh));
   %Fst(dh)=pi*a*Kb*T/s^3*((8*Ii/5)*((2*Ii/h(dh))^(5/4)-1)+(8*Ii/7)*((h(dh)/(2*Ii))^(7/4)-1));
   %Vtot(dh)=(Vr(dh)-Va(dh))/(Kb*T);
   if h(dh)>=2*Ii
      Vo(dh)=0; Ve(dh)=0; 
   elseif h(dh)>=Ii && h(dh)<2*Ii
      Vo(dh)=4*pi*a/(0.01*2.99e-27)*phip^2*(1/2-qsi)*(Ii-h(dh)/2)^2*(Kb*T);
     % Vo(dh)=4*pi*(a*1e9)/0.03*phip^2*(1/2-qsi)*(Ii*1e9-h(dh)*1e9/2)^2*Kb*T;
      Ve(dh)=0;
   else
      Vo(dh)=4*pi*a/(0.01*2.99e-27)*phip^2*(1/2-qsi)*Ii^2*(h(dh)/(2*Ii)-1/4-log(h(dh)/Ii))*(Kb*T); 
     % Vo(dh)=4*pi*(a*1e9)/0.03*phip^2*(1/2-qsi)*(Ii*1e9)^2*(h(dh)/(2*Ii)-1/4-log(h(dh)/Ii))*(Kb*T);
      Ve(dh)=2*pi*(a)*6.023e23/(Mw)*phip*(Ii)^2*rho_lig*1e6*...
          (h(dh)/Ii*log(h(dh)/Ii*((3-h(dh)/Ii)/2)^2)-6*log((3-(h(dh)/Ii))/2)+3*(1-(h(dh)/Ii)))*(Kb*T);
      
      
   end
   
   Vtot(dh)=(Vedl(dh)-Vvdw(dh)+0*Fst(dh)+Vo(dh)+Ve(dh))/(Kb*T);
   
dh=dh+1;    
end


dh=numel(h);
while dh>2
    if dh<numel(h)
    Vst(dh)=Vst(dh+1)+Fst(dh)*(h(2)-h(1))/(Kb*T);    
    end
dh=dh-1;
end
Vst(1)=Vst(2);
 %plot(h,Vtot);
 %[Vmax, hmax]=findpeaks(Vtot);


%txt = '\leftarrow Energy barrier';
%text(h(hmax),Vmax,txt,'FontSize',14);
%[Vmin, hmin]=findpeaks(-1*Vtot);
%if numel(Vmin)==2
%Vmin=Vmin; hmin=hmin;
%ylim([-10 abs(max(Vmax))*2]); set(gcf,'color','w');
%txt = '\leftarrow Secondary minimum';
%text(h(hmin),-Vmin,txt,'FontSize',14);
%elseif numel(Vmax)>0
%ylim([-10 abs(max(Vmax))*2]); set(gcf,'color','w');    
%else
%ylim([-10 0]); set(gcf,'color','w');     
%end

%txt = '\leftarrow Primary minimum';
%text(0,-10,txt,'FontSize',14);
% VMAX=abs(Vmax-Vmin);
Wm=1/(2*Cinv*a)*exp(0);
dh=1;
while dh<=numel(h)
    u=h(dh)/a;
bh(dh)=(6*(u)^2+13*(u)+2)/(6*(u)^2+4*(u));
   if dh>1
   %int_Vtot(dh)=int_Vtot(dh-1)+bh(dh)*exp((Vr(dh)-Va(dh))/(Kb*T))/(2+u)^2*(h(2)-h(1))/a;
   int_Vtot(dh)=int_Vtot(dh-1)+bh(dh)*exp((Vedl(dh)+Vo(dh)+Ve(dh)-Vvdw(dh))/(Kb*T))/(2+u)^2*(h(2)-h(1))/a;
   int_Vvdw(dh)=int_Vvdw(dh-1)+bh(dh)*exp(-1*Vvdw(dh)/(Kb*T))/(2+u)^2*(h(2)-h(1))/a;
   end
dh=dh+1;
end
% if Vmax>0
% W=1/(2*Cinv*a)*exp(Vmax);
% alpha=exp(-Vmax);
% else
% alpha=1;
% end
%f=sum(Vo);
f=1/(max(int_Vtot)/max(int_Vvdw));
%f=abs(1/(max(int_Vtot)/max(int_Vvdw))-fthe)/fthe;
%f=Vtot