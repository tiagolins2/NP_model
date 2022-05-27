
%tt=0;    
%dt=t_end/1000;
%Results_level3=[];
iii=1;

%while tt<=t_end 
Vair=9.42e14;
Vwater=1.63e12;
Vsoil=7.53e9;
Vsed=1.9e8+((f_sedtotal(end))*4/3*pi*(arr_sizes(1)/2)^3+sum(sum((4/3*pi*(arr_sizes_ss/2).^(3)).*nssij_sed)))*Vwater;
Vaer=Vair*2e-11;
%Vsus=3.4e5+f_hete(end)*4/3*pi*(arr_sizes(1)/2)^3*Vwater;
Vsus=(sum((4/3*pi*(arr_sizes_ss/2).^(3)).*nss)+sum(sum((4/3*pi*(arr_sizes_ss/2).^(3)).*nssij')))*Vwater;
Vfish=1.63e5;
%avg_size=sum((y(last_d,1:numel(arr_sizes)))'.*arr_sizes)/sum(y(last_d,1:numel(arr_sizes)));
%num_of_part=sum(y(last_d,1:numel(arr_sizes)));
%Vnp=((f_homo(end)+f_single(end))*4/3*pi*(arr_sizes(1)/2)^3+sum(sum((4/3*pi*(arr_sizes/2).^(3)).*nssij')))*Vwater;
%Vnp=((num_of_part)*4/3*pi*(avg_size/2)^3+sum(sum((4/3*pi*(arr_sizes/2).^(3)).*nssij')))*Vwater;
%Anp=((f_homo(end)+f_single(end))*4*pi*(arr_sizes(1)/2)^2+sum(sum((4*pi*(arr_sizes/2).^(2)).*nssij')))*Vwater;
%Anp=((num_of_part)*4*pi*(avg_size/2)^2+sum(sum((4*pi*(arr_sizes/2).^(2)).*nssij')))*Vwater;
Awater=1.9e10;
Asoil=7.53e10;
Ap=4*pi*(arr_sizes(1)/2)^2;
Vp=4/3*pi*(arr_sizes(1)/2)^3;
Knp=3.04;
Kpw=10^Knp;
%del_w=arr_sizes(1)/2;
del_w=50e-6;
%del_p=(sum(n.*arr_sizes)/sum(n))/2;
del_p=avg_size/2;
D_x_w=1e-9;
mp=7.41e-16*Kpw/del_p;
f_np_v=Vnp/Vwater;
k_u=Anp/Vnp*(del_w/D_x_w*(avg_size/2)/(del_w+avg_size/2)+1/mp)^-1;%*f_np_v;
k_r=k_u/Kpw;
%Vwnp=4/3*pi*(((del_w*2+arr_sizes(1))/2)^3-(avg_size/2)^3)*(f_homo(end)+f_single(end))*Vwater+sum(sum((4/3*pi*(arr_sizes_ss/2).^(3)).*nssij'))*Vwater;
%if Vwnp>Vwater Vwnp=Vwater; end

V_arr=[Vwater;Vair;Vsoil;Vsed;Vaer;Vsus;Vfish;Vnp];

yfish=0.05;
ysus=0.2;
ysoil=0.02;
ysed=0.0359;

Pwater=rho_w;
Pair=1.185;
Psoil=2400;
Psed=2400;
Paer=2000;
Psus=2400;
Pfish=1000;
Pnp=rho_p;

% Molar_mass=252.13; %228.29 252.31 228.29 252.31 252.31 276.33
% %temp=9.48;  
% Melting_pt=175; % 254 163 158 175 217 278
% VP=5.74E-5;     %1.2e-7  6.67e-5  2.8e-5  5.74E-5,7.3e-7  1.3e-8  2.55e-5
% CW=0.0038;   %2e-3    1.5e-3   9.4e-3  3.8e-3   8e-4   0.079
%        
% H=0.046445;  %0.5299297, 0.06657053, 1.2159, 0.04630553, 0.0591738, 0.03353858
%              %CHR BBF BAA BAP BKF BGP
% logKow=6.04; %5.8 5.8 5.9 6.0 6.2 6.9
% Koc=0.41*10^logKow; 
% Koc=4.04e5;
% F_ratio=exp(6.79*(1-(Melting_pt+273.15)/(T)));
% SLVP=VP/F_ratio;
% M=222.26;
%            % CHR BBF  BAA BAP  BKF  BGP
% Knp=7.95;  %7.62 7.65 7.3 7.95 7.61 8.04
Z1=1/H;
Z2=1/(8.314*(T));
Z3=Z1*Psoil*ysoil*Koc/1000;
Z4=Z1*Psed*ysed*Koc/1000;
Zaer=6000000*Z2/SLVP;
Zsus=Z1*Psus*ysus*Koc/1000;
Zfish=Z1*yfish*10^logKow;
Znp=Z1*10^Knp*densityNP_a/1000;   %
Zwater=(1-Vsus/Vwater-Vfish/Vwater-Vnp/Vwater)*Z1+Vsus/Vwater*Zsus+Vfish/Vwater*Zfish+Vnp/Vwater*Znp;
Zair=(1-Vaer/Vair)*Z2+Vaer/Vair*Zaer;
Zsoil=0.2*Z2+0.3*Z1+0.5*Z3;
Zsed=0.8*Z1+0.2*Z4+Znp*(f_sedtotal(end))*4/3*pi*(arr_sizes(1)/2)^3;
Z_arr=[Zwater;Zair;Zsoil;Zsed;Zaer;Zsus;Zfish;Znp];

%                     % chrysene  BBF    BAA      BAP       BKF         BGP
% Cwi=10.55;  %pg/L    151.75     44.97  227.88   10.55     34.15       912.8
% Cwi_mol=Cwi*1e-9/Molar_mass; 
% Cai=25;     %pg/m3   220        52     170      25        43          45    
% Cai_mol=Cai*1e-12/Molar_mass;
% Csedi=0.025; %ug/g              0.049  0.017    0.025     0.049       0.05?
% Csedi_mol=Csedi*1e-6*1000*Psed/Molar_mass;
Tair=23.4*2;
Twater=6.89e4;
Tsed=256000;
MTC_21=1;
MTC_22=0.01;
MTC_23=9.59e-5;
MTC_24=0.72;
MTC_25=0.02;
MTC_26=1e-5;
MTC_27=5;
MTC_28=0.0002;
MTC_29=4.557e-8;
MTC_30=1.9531e-8;
MTC_31=0.00005;
MTC_32=0.00000001;
D125=Awater*MTC_21*Z2;
D126=Awater*MTC_22*Z1;
D127=1/((1/D125)+(1/D126));
D129=Awater*MTC_23*Z1;
D130=Asoil*MTC_23*Z1;
D132=Awater*Vaer/Vair*MTC_24*Zaer;
D133=Awater*Vaer/Vair*MTC_23*Zaer*200000;
D134=D132+D133;
D136=Asoil*Vaer/Vair*MTC_24*Zaer;
D137=Asoil*Vaer/Vair*MTC_23*Zaer*200000;
D138=D136+D137;
D125h=Asoil*MTC_31*Z1;
D126h=Asoil*MTC_32*Z3;
D128h=Asoil*MTC_25*Z2;
D129h=Asoil*MTC_26*Z1;
D130h=Asoil*MTC_27*Z2;
D131h=1/((1/D128h)+(1/(D129h+D130h)));
D133h=Awater*MTC_28*Z1;
D135h=Awater*MTC_29*Zsus*0;
D136h=Awater*MTC_30*Z4*0;

D21T=D127+D134+D129;
D23T=D131h+D130+D138;
D12T=D127;
D14T=D133h+D135h;
D32T=D131h;
D31T=D125h+D126h;
D41T=D133h+D136h;
% D21=MTC21*Awater*Zair;
% D12=MTC12*Awater*Zwater;
% D23r=MTC21r*Asoil*Zwater;
% D21r=MTC21r*Awater*Zwater;
% D23=MTC23*Asoil*Zair;
% D23b=MTC23b*Asoil*Zair;
% D31=MTC31*Asoil*Zwater;
% D41d=MTC41*Awater*Zwater;
% D14=MTC14*Awater*Zsus;
% D41=MTC41*Awater*Zsed;
% D31r=MTC31r*Asoil*Zwater;
% D36r=MTC36r*Asoil*Zsoil;
% D21d=MTC21d*Awater*Zaer;
% D51w=Awater*Zaer*MTC21r*Vaer/Vair;
% D23d=MTC21d*Asoil*Zaer;
% D53w=Asoil*Zaer*MTC21r*Vaer/Vair;
% D21T=(D21+D21d+D51w+D21r);
% D12T=D12;
% D23T=(D23r+D23d+D53w+1/((1/D23)+1/(D23b+D31)));
% D14T=D14+D41d;
% D41T=D14+D41;
% D32T=1/((1/D23)+1/(D23b+D31));
% D31=D31r+D31
% D21T=1.183e9;
% D23T=4.665e9;
% D12T=7.634e6;
% D14T=8.164e7;
% D32T=6.051e5;
% D31T=8.103e7;
% D41T=8.164e7;
D2R=7.126e9;
D1R=1.049e10;
D3R=7.128e10;
D4R=2.762e8;
%D1A=5.605e8;  
D1A=Vwater*(1-Vsus/Vwater-Vfish/Vwater-Vnp/Vwater)*Z1/Twater+...
    Zsus*Vsus/Twater+Zfish*Vfish/Twater+Znp*Vnp/Twater;
D2A=Z2*Vair/Tair+Zaer*Vaer/Tair;
D3A=0;
D4A=1.235E8;
Ew=3900*0.02/0.98*1e3/(365*24*Molar_mass)*0.1;
Ew=Ew+Vwater/Twater*Cwi_mol;
Ea=3900e3/(365*24*Molar_mass)*0.1;
Ea=Ea+Vair/Tair*Cai_mol;
%kr=2.13e-5;
%ku=kr*Kpw;

% m=[(D12T+D14T+D1R+D1A+k_u) -D21T -D31T -D41T -k_r*f_np_v;...
%     -D12T (D21T+D23T+D2R+D2A) -D32T 0 0; 0 -D23T (D32T+D31T+D3R+D3A) 0 0;...
%     -D14T 0 0 (D41T+D4R+D4A) 0;-k_u 0 0 0 k_r*f_np_v];
m=[(D12T+D14T+D1R+D1A+k_u*Vnp) -D21T -D31T -D41T -k_r*Vnp;...
    -D12T (D21T+D23T+D2R+D2A) -D32T 0 0; 0 -D23T (D32T+D31T+D3R+D3A) 0 0;...
    -D14T 0 0 (D41T+D4R+D4A) 0;-k_u*Vnp 0 0 0 k_r*Vnp];
I=[Ew; Ea; 0; 0; 0];
%Fugacity=m\I;
% Fugacity=([(D12T+D14T+D1R+D1A) -D21T -D31T -D41T;...
%     -D12T (D21T+D23T+D2R+D2A) -D32T 0; 0 -D23T (D32T+D31T+D3R+D3A) 0;...
%     -D14T 0 0 (D41T+D4R+D4A)])\I;

% Fugacity(1)=(15E-9*1000/Molar_mass)/Z1;
% Fugacity(2)=(0.05E-9/Molar_mass)/Z2;
% Fugacity(3)=(782e-6/(Molar_mass*Psed))/Z3;
% Fugacity(4)=(782e-6/(Molar_mass*Psed))/Z4;

%fg=Fugacity';
% dfdt=([(D12T+D14T+D1R+D1A) -D21T -D31T -D41T;...
%     -D12T (D21T+D23T+D2R+D2A) -D32T 0; 0 -D23T (D32T+D31T+D3R+D3A) 0;...
%     -D14T 0 0 (D41T+D4R+D4A)]*Fugacity'-[Ew; Ea; 0; 0])/(V_arr(1:4).*Z_arr(1:4));

VZ=[V_arr(1:4);V_arr(8)].*[Z_arr(1:4);Z_arr(8)];
%t_end=1000;
%tt=2;
K_in=(1/5*Vfish*(Vnp/Vwater))/(12); %m3np/h
K_out=log(2)/48; %1/h
%fg(5)=K_in/(K_out*Vfish);
%dCdt=(K_in-K_out*C_fish*Vfish)/Vfish;

%options = optimset('Display','off');
%[tu,yy]=ode23(@(tu,fg)level(fg,m,I,VZ,K_in,K_out,Vfish),[tt/(3600) (tt+dt)/(3600)]',fg,options);

%Fugacity=yy(numel(yy)/5,:)
VZ_arr=V_arr.*Z_arr;
%fugacity=M/sum(VZ_arr);
% C1=Fugacity(1)*Z1;
% C2=Fugacity(2)*Z2;
% C3=Fugacity(3)*Z3;
% C4=Fugacity(4)*Z4;
% C5=Fugacity(2)*Z_arr(5);
% C6=Fugacity(1)*Z_arr(6);
% C7=Fugacity(1)*Z_arr(7);
% C8=Fugacity(1)*Z_arr(8);
% C_arr=[C1;C2;C3;C4;C5;C6;C7;C8];
% m_arr=C_arr.*V_arr;
%C_arr=Z_arr*fugacity;
%m_arr=VZ_arr*fugacity;

%Results_level3=[Results_level3; Vnp/Vwater m_arr' C_arr' Fugacity(5) ]; 
%tt=tt+dt;
iii=iii+1;
%end