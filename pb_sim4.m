%function f=pb_sum2(NPsize_a,NCsize_a,NCconc_a,shear_a,densityNP_a,densityNC_a,densityNOM_a,salinity_a,fractal_d_a,Hamaker_A_a,Hamaker_Ass_a,depth_a,Rmax_a, OM_size_a, OM_packing_a, Temperature_a)
ccc=0;
ff=[];
ff_het=[];
ff_hom=[];
ff_single=[];
ff_sedtotal=[];
ff_np_in_fish=[];
ff_fish=[];
ff_water=[];
ff_np=[];
ff_sed=[];
ff_fish_np=[];
ff_total_np_system=[];
ff_air=[];
ff_sediment=[];
ff_sus_sediment=[]; ff_avg_size=[];
ff_Knp=[];
%
arr_K=zeros(71,5,5);
arr_Kss=zeros(71,5,5); Fugacity=[0 0 0 0];
while ccc<=50

size_p=100*1e-9;
%arr_sizes=size_p*([1 2 3 4 5 6 7 8 9 10]).^(1/3);
%arr_sizes=size_p*([1 2 4 8 16 32 64]).^(1/3);
arr_sizes=2.^(linspace(0,9,10));
arr_sizes_ss=[1 1 2 2 4 4 8 8 16 16];
%NCsize_a=100000;
%NPsize_a=103.5;
%arr_sizes=linspace(1,400,401);
%arr_sizes=arr_sizes.^(1/3);
arr_sizes=arr_sizes*size_p;
arr_sizes=arr_sizes';
arr_sizes_ss=arr_sizes_ss*size_p;
arr_sizes_ss=arr_sizes_ss';
arr_sizes_ss=arr_sizes_ss*NCsize_a/100;
arr_sizes=arr_sizes*NPsize_a/100*1e-2*1.3183^ccc;

%depth_a=0.1;
n=zeros(numel(arr_sizes),1);
%NCconc_a=100*0;
%densityNC_a=100;
nss=ones(numel(arr_sizes),1)*(0.001*NCconc_a/(sum(4/3*pi()*(arr_sizes_ss/2).^3)*densityNC_a));%*500^3/NCsize_a^3;%*(1e-3*1.5^ccc);
nssij=zeros(numel(arr_sizes),numel(arr_sizes));

n_sed=zeros(numel(arr_sizes),1);
nss_sed=ones(numel(arr_sizes),1)*1.54e9*10000;
nssij_sed=zeros(numel(arr_sizes),numel(arr_sizes))*1e2;
n_np_in_fish=zeros(numel(arr_sizes),1);
%inital_c=(5.73e13);
%inital_c=1.00605e16;
%inital_c=1.82e12*1e-3*1.5^ccc;
%concentrationNP_a=10000;
inital_c=concentrationNP_a*1e-6/(densityNP_a*4/3*pi*(arr_sizes(1)/2)^3); 
%inital_c=inital_c*(100*1e-9)^3/arr_sizes(1)^3;
n(1,1)=inital_c;
n0=n;
t=0;
K=zeros(numel(n),numel(n));
Kss=zeros(numel(n),numel(nss));
hydro_size=[];
v_dndt=zeros(numel(arr_sizes),1);
kb=1.380648813131313131313131e-23;
g=9.81;
T=Temperature_a;
mu=-9.80215e-9*T^3+9.16232E-6*T^2-2.86927E-3*T+3.01685E-1;
s=0;
I=0; %I=100+ccc*25; 
%salinity_a=800;
I=salinity_a;%+ccc*20;
s=I/1000*58.44;
df=fractal_d_a;
mu=mu*(1.9609e-3*s+9.9276e-1);
rho_p=densityNP_a;
rho_ps=densityNC_a;
rho_w=1;
%shear_perc=0.1;

rho_w=rho_w*(0.7223*s+4.497e-5*T^3-4.518E-2*T^2+1.469E1*T-5.576E2);
rhow=rho_w;
t_end=12*30*24*60*60;
t_end=60*60*24*30;% t_end=60*60*4*ccc; t_end=600*1;
G=Shear_a; 
%G=shear_a;
Gc=G*mu/0.5; tau_v=Gc^2*rhow; shear_perc=tau_v/0.3-1; %G=0;
shear_perc=0e-9;
if (shear_perc<0) shear_perc=0; end 
%shear_perc=10;
%t_end=10;
%tau_v=(shear_perc+1)*0.5;
%Gc=sqrt(tau_v/rhow);
%G=0.5*Gc/mu*0;
%biofilm/DOM parameters:
phi_pa=OM_packing_a; L=OM_size_a; Mw=500000;
rho_pa=densityNOM_a*1000;
rho_L=rho_pa*phi_pa+rhow*(1-phi_pa);
V_L=4/3*pi*((size_p/2+L)^3-(size_p/2)^3);
rho_pL=(rho_p*4/3*pi*(size_p/2)^3+rho_L*V_L)/(4/3*pi*((size_p/2+L)^3));
arr_sizes=arr_sizes*(size_p+L*2)/size_p; %adjusted size to take into account DOM layer
rho_p=rho_pL; %adjusted density based on packing of organic ligands.

%adding particle dynamic effect for NOM/DOM
if active_stoch==1
n_dom=ccc*10;
Xi=50;
X=Xi; %lower is faster X_2=0;
X_4=0;X_8=0;X_16=0;X_32=0;X_64=0;Xdom=X*n_dom; Xdomi=Xdom; nn_domfrac=0;
alpha_np_np=alpha; alpha_dom_dom=0; alpha_dom_np=1;
Rdom=1.0e-9/2; C=inital_c/(1.82e15); %mg/L
DOM_size=2*Rdom; DOM_NumPerNP=nn_domfrac; rho_dom=1330;
end
%alpha=0.01;
%dom_size_arr=[1E-10	0.000000001	0.000000002	0.000000003	0.000000004	0.000000005	0.000000006	0.000000007	0.000000008	0.000000009	0.00000001	0.000000011	0.000000012	0.000000013	0.000000014	0.000000015	0.000000016	0.000000017	0.000000018	0.000000019	4.00E-08];
%prob_dd_arr=[0	0	0.000189448	0.000686292	0.001788269	0.003884044	0.006490089	0.008832565	0.014289669	0.018640763	0.027661902	0.038332372	0.05355581	0.068674256	0.093282712	0.115289542	0.150434925	0.172593833	0.210936735	0.238910177	0.853054402];
%prob_nn_arr=[0.995175608	0.946558808	0.882195068	0.838324002	0.770958512	0.710179725	0.663111735	0.651270772	0.593003944	0.558573041	0.526670314	0.495912344	0.461010563	0.409575001	0.383354981	0.328070927	0.299458395	0.301086921	0.271036575	0.249886341	0.024617567];
%prob_dd_arr=[0	0.001979967	0.019730839	0.06031746	0.124459343	0.226960994	0.36112592	0.500319285	0.658619589	0.776726886	0.862046976	0.931038147	0.958018918	0.989159316	0.99479222	0.997662807	0.998724761	0.999787483	0.999468594	0.999681258	0.999747666];
%prob_nn_arr=[0.998963596	0.490798975	0.326528285	0.147016968	0.104346886	0.044470639	0.019191811	0.011813538	0.006699989	0.001487779	0.001062812	0.001700138	0.000637687	0.000106281	0	0.000106236	0.00010627	0.000212517	0.000318844	0.000212495	0];
%prob_dd_arr=[0	0.161569533	0.761520737	0.984154422	0.999711816	0.999711816	0.998559908	1	0.999711899	0.999711732	0.999423631	1	0.999423465	0.998847594	0.999423631	1	1	0.999711732	1	1	0.999747666];
%prob_nn_arr=[0.997734994	0.005481823	0.000576037	0.000576203	0.000288184	0.000288184	0.000864055	0	0.000288101	0	0	0	0	0.000288101	0	0	0	0	0	0	0];



% DOM_size=dom_size_arr(ccc);
% DOM_NumPerNP=1000;
% prob_dom_dom=prob_dd_arr(ccc);
% prob_np_np=prob_nn_arr(ccc);
% prob_np_dom=1-(prob_dom_dom+prob_np_np);
% alpha=alpha*prob_np_np+prob_np_dom;
% rho_dom=1080;
if active_stoch==1
DOM_size=2*Rdom;  
Vt=4/3*pi*(size_p/2+DOM_size)^3;
Vnp=4/3*pi*(size_p/2)^3; 
V_dom=DOM_NumPerNP*4/3*pi*(DOM_size/2)^3; Vw=Vt-(Vnp+V_dom);
rho_domNP=(rhow*Vw+rho_p*Vnp+rho_dom*V_dom)/Vt;
arr_sizes=arr_sizes*(size_p+DOM_size*2*(DOM_NumPerNP>0))/size_p;
rho_p=rho_domNP;
end

%alpha=DVLO(I);
%alpha_ss=DVLO(I);



%alpha=0.001;

alpha=pb_xDLVO(I, arr_sizes(1)/2, T, 0.485, phi_pa, L, rho_pa,Mw,(sqrt(Hamaker_A_a)-sqrt(3.7e-20))^2);
%alpha_ss=alpha_ss_adj;

alpha_ss=pb_het_DLVO(I,T,0.485, phi_pa, L, rho_pa, Mw,Hamaker_A_a,Hamaker_Ass_a,arr_sizes(1)/2,arr_sizes_ss(1)/2);
if alpha>0.004
   df=0.17*log10(1/alpha)+1.67;
else
   df=2.1; 
end
%df=2.5;
%alpha=10^-3.75;
%alpha=1;
%alpha_ss=alpha;
if isnan(alpha) alpha=1; end
if isnan(alpha_ss) alpha_ss=1; end

%alpha=0.1;
%alpha_ss=0.0;
alpha_r=alpha;
alpha_ss_r=alpha_ss;
%alpha_r=alpha*(1-(0.5-rand()));
%alpha_ss_r=alpha_ss*(1-(0.5-rand()));


t=0;
i=1;
j=1;
rhop=zeros(numel(n),1);
rhops=zeros(numel(n),1);
a1=arr_sizes(1)/2;
dt=10;
while i<=numel(n)
    ai=arr_sizes(i)/2;
    rhop(i)=((ai/a1)^df*4/3*pi*a1^3*rho_p+rho_w*(4/3*pi*ai^3-(ai/a1)^df*4/3*pi*a1^3))/(4/3*pi*ai^3);
    rhops(i)=((ai/a1)^df*4/3*pi*(a1*5)^3*rho_ps+rho_w*(4/3*pi*(ai*5)^3-(ai/a1)^df*4/3*pi*(a1*5)^3))/(4/3*pi*(ai*5)^3);
    i=i+1;
end

i=1;
while i<=numel(n)
    j=1;
    while j<=numel(n)
    %assuming no shear rate:
    ai=arr_sizes(i)/2*i^(1/df);
    aj=arr_sizes(j)/2*j^(1/df);
    %ai=arr_sizes(i)/2;
    %aj=arr_sizes(j)/2;
    K(i,j)=((2*kb*T)/(3*mu)*(ai+aj)^2/(ai*aj)+4/3*G*(ai+aj)^3+...
        2*pi*g/(9*mu)*(rhop(j)-rho_w)*(ai+aj)^3*(ai-aj));
    j=j+1;
    end
    i=i+1;
end
    
i=1;
j=1;
while i<=numel(nss)
    j=1;
    while j<=numel(n)
        ai=arr_sizes_ss(i)/2;
        aj=arr_sizes(j)/2;
        vsj=2*aj^2*(rhop(j)-rhow)*g/(9*mu);
        vsi=2*ai^2*(rho_ps-rhow)*g/(9*mu);        
        Kss(i,j)=((2*kb*T)/(3*mu)*(ai+aj)^2/(ai*aj)+4/3*G*(ai+aj)^3+...
        pi*(aj+ai)^2*abs(vsj-vsi));
        j=j+1; 
    end
    i=i+1;
end


options = optimset('Display','off');
options.OptimalityTolerance=1e-22;
options.FunctionTolerance=1e-11;
ii=1;
%while t<t_end 
%i=1;

% j=1;
% while j<numel(n) 
%     i=1;
%     s1=0;
%     s2=0;
%     while i<numel(n)
% %         s2=s2+alpha*K(i,j)*n(i);
% %         if j==1
% %             s1=0;
% %         end
% %         if j>1 && i<=j-1
% %             s1=s1+alpha*K(i,j-1)*n(i)*n(j-1);
% %         end
%         dndt(K(i,j), n(i), K(i,j-1),n(j-1),alpha,j>1,i<=j-1)
%         i=i+1;
%     end
%     
%     v_dndt(j)=1/2*s1-n(j)*s2; 
%       
%    j=j+1; 
% end

tt=0;    
dt=t_end/100;

fnp=[];
 f_hete=[0];
 f_single=[inital_c];
 f_homo=[0];
 f_fish_np=[0];
 f_avg_size=[arr_sizes(1)];
 f_total_np_system=[inital_c];
 f_sedtotal=[0];
 avg_size=sum(n.*arr_sizes)/sum(n);
 Vnp=((sum(n.*arr_sizes)/arr_sizes(1))*4/3*pi*(arr_sizes(1)/2)^3+sum(sum((4/3*pi*(arr_sizes/2).^(3)).*nssij')))*Vwater;
 Anp=((sum(n.*arr_sizes)/arr_sizes(1))*4*pi*(arr_sizes(1)/2)^2+sum(sum((4*pi*(arr_sizes/2).^(2)).*nssij')))*Vwater;

 level4;
 
 
 Fugacity(1)=Cwi_mol/Z1;
 Fugacity(2)=Cai_mol/Z2;
 Fugacity(3)=(Csedi_mol*4)/Z3; 
 Fugacity(4)=Csedi_mol/Z4;
 Fugacity(5)=0;
  
fg=Fugacity; %fg(5)=K_in/(K_out*Vfish); %fg(5)=0;
 fg=fg';
 y_val=[n;nss;nssij(:);fg;n_sed;nss_sed;nssij_sed(:)]';
 Results_level3=[];


 tt=0;
avg_Knp=[]; 
Knp_arr_process=[];
Knp_process=[];
while tt<=t_end 
    n_prev=n;
    
    n0 = n_prev;
    
    
    
    % add level4 to n0:
    n0=[n;nss;(nssij(:))];
    avg_size=sum(n.*arr_sizes)/sum(n);
 Vnp=((sum(n.*arr_sizes)/arr_sizes(1))*4/3*pi*(arr_sizes(1)/2)^3+sum(sum((4/3*pi*(arr_sizes/2).^(3)).*nssij')))*Vwater;
 Anp=((sum(n.*arr_sizes)/arr_sizes(1))*4*pi*(arr_sizes(1)/2)^2+sum(sum((4*pi*(arr_sizes/2).^(2)).*nssij')))*Vwater;

    level4;
    n0=[n0;fg];
    n0=[n0;n_sed;nss_sed;(nssij_sed(:))];
    
    ntotal=n0;
    %options = odeset('RelTol',1e-30,'Stats','on');
    nu_sizes=numel(arr_sizes);
    %[t,ysimple]=ode15s(@(t,n)hetero_odesimple(n,nss,nssij,ntotal,n_prev,K,Kss,dt,alpha,alpha_ss,arr_sizes,t,rhow,rhop,rho_ps,mu,numel(arr_sizes)),[0 t_end]',[inital_c 0 0 0 0]);
    %[t,y]=ode15s(@(t,ntotal)hetero_ode(ntotal,n_prev,K,Kss,dt,alpha,alpha_ss,arr_sizes,t,rhow,rhop,rho_ps,mu,numel(arr_sizes)),[tt tt+dt]',ntotal);
    %rand_size
    %(1-(0.5-rand())*0.2);
    %rand_density
    %rand_alpha

    %rand_concentration
    %rand
     
%      if (tt>t_end/10)
%          shear_perc=1e-11;
%      else
%          shear_perc=0;
%      end

   if active_stoch==1 
   stoch_pbsim; %alters alpha_r  
   
   %correct size and density of particles
   DOM_NumPerNP=nn_domfrac; rho_dom=1330;
Vt=4/3*pi*(size_p/2+DOM_size)^3;
Vnp=4/3*pi*(size_p/2)^3; 
V_dom=DOM_NumPerNP*4/3*pi*(DOM_size/2)^3; Vw=Vt-(Vnp+V_dom);
rho_domNP=(rhow*Vw+rho_p*Vnp+rho_dom*V_dom)/Vt;
rho_p=rho_domNP;
arr_sizes=2.^(linspace(0,4,5))*(size_p+2*DOM_size*(DOM_NumPerNP>0)); 
arr_sizes=arr_sizes';
end
   %
   %estimate rates at each time-step:
   rates_nps;
Knp_arr_process=[Knp_arr_process; rates_of_np_arr];
Knp_process=[Knp_process; rates_of_np];

%run the solver:
    [t,y]=ode15s(@(t,ntotal)pb_ode4(ntotal,n_prev,K,Kss,dt,alpha_r,alpha_ss_r,arr_sizes,arr_sizes_ss,t,rhow,rhop,rho_ps,mu,numel(arr_sizes),m/3600,I/3600,VZ,[Z1 0 0 Z4],K_in/3600,K_out/3600,Vfish,shear_perc, inital_c,Rmax_a*0, depth_a,resus_rate_a,df),[tt tt+dt]',ntotal);
   % [t,y]=ode15s(@(t,ntotal)homo_ode_fugacity(ntotal,n_prev,K,Kss,dt,alpha_r,alpha_ss_r,arr_sizes,arr_sizes_ss,t,rhow,rhop,rho_ps,mu,numel(arr_sizes),m/3600,I/3600,VZ,K_in/3600,K_out/3600,Vfish,shear_perc, inital_c),[tt tt+dt]',ntotal); 
    [last_d,mmm]=size(y);
    
    %sum((y(last_d,:))'.*arr_sizes)/sum(y(last_d,:));
    %plot(arr_sizes(1:5),y(last_d,1:5));
    jj=1;
    arr_results_lv3=[];
    
   % while jj<last_d
        
    fnp=[fnp; sum((y(last_d,1:numel(arr_sizes)))'.*arr_sizes)/(arr_sizes(1)*inital_c)*inital_c*4/3*pi*(arr_sizes(1)/2)^3];
    
    
    
    avg_size=sum((y(last_d,1:numel(arr_sizes)))'.*arr_sizes)/sum(y(last_d,1:numel(arr_sizes)));
    aggr_sum=0;
    for jjn=1:numel(arr_sizes)
    aggr_sum=aggr_sum+sum((y(last_d,numel(arr_sizes)+numel(arr_sizes)*jjn+1:numel(arr_sizes)+numel(arr_sizes)*jjn+numel(arr_sizes)))'*arr_sizes(jjn));
        %aggr_sum=sum((y(last_d,21:20+numel(arr_sizes)))'*arr_sizes(1))+sum((y(last_d,31:30+numel(arr_sizes)))'*arr_sizes(2))+sum((y(last_d,41:40+numel(arr_sizes)))'*arr_sizes(3))+sum((y(last_d,51:50+numel(arr_sizes)))'*arr_sizes(4))+sum((y(last_d,61:60+numel(arr_sizes)))'*arr_sizes(5))+...
        %sum((y(last_d,71:70+numel(arr_sizes)))'*arr_sizes(6))+sum((y(last_d,81:80+numel(arr_sizes)))'*arr_sizes(7))+sum((y(last_d,91:90+numel(arr_sizes)))'*arr_sizes(8))+sum((y(last_d,101:100+numel(arr_sizes)))'*arr_sizes(9))+sum((y(last_d,111:110+numel(arr_sizes)))'*arr_sizes(10));
    end
    aggr_sum=aggr_sum/(arr_sizes(1)*inital_c);
    single_part_sum=y(last_d,1)*arr_sizes(1)/(arr_sizes(1)*inital_c);
    homo_sum=sum(y(last_d,2:nu_sizes).*(arr_sizes(2:nu_sizes))')/(arr_sizes(1)*inital_c);
    %sed_particles=sum((y(last_d,11:10+numel(arr_sizes)))'*arr_sizes(1))+sum((y(last_d,16:15+numel(arr_sizes)))'*arr_sizes(2))+sum((y(last_d,21:20+numel(arr_sizes)))'*arr_sizes(3))+sum((y(last_d,26:25+numel(arr_sizes)))'*arr_sizes(4))+sum((y(last_d,31:30+numel(arr_sizes)))'*arr_sizes(5));
    
    %f_total_np_system=[f_total_np_system; aggr_sum*inital_c+single_part_sum*inital_c+homo_sum*inital_c];
    f_hete=[f_hete; aggr_sum*inital_c];
    f_single=[f_single; single_part_sum*inital_c];
    f_homo=[f_homo; homo_sum*inital_c];
    %f_fish_np=[f_fish_np; y(last_d,76)];
    f_avg_size=[f_avg_size; avg_size];
    
    
    % level_3
    
   % arr_results_lv3=[arr_results_lv3; [t(jj) fnp Results_level3]];
    jj=jj+1;
   % end
    n=y(last_d,1:nu_sizes);
    n=n';
    nss=y(last_d,nu_sizes+1:2*nu_sizes);%(nu_sizes+1:2*nu_sizes);
    nss=nss';
    nssij=y(last_d,2*nu_sizes+1:2*nu_sizes+nu_sizes*nu_sizes);
    nssij=reshape(nssij,nu_sizes,nu_sizes);
    fg=y(last_d,2*nu_sizes+nu_sizes*nu_sizes+1:2*nu_sizes+nu_sizes*nu_sizes+5);
    fg=fg';
    
    n_sed=y(last_d,2*nu_sizes+nu_sizes*nu_sizes+6:2*nu_sizes+nu_sizes*nu_sizes+5+nu_sizes);
    n_sed=n_sed';
    nss_sed=y(last_d,2*nu_sizes+nu_sizes*nu_sizes+6+nu_sizes:2*nu_sizes+nu_sizes*nu_sizes+5+2*nu_sizes);
    nss_sed=nss_sed';
    nssij_sed=y(last_d,2*nu_sizes+nu_sizes*nu_sizes+6+2*nu_sizes:2*nu_sizes+nu_sizes*nu_sizes+5+2*nu_sizes+nu_sizes*nu_sizes);
   
    nssij_sed=reshape(nssij_sed,nu_sizes,nu_sizes);
    %n_np_in_fish=y(last_d,76:80); n_np_in_fish=n_np_in_fish';
    y_val=[y_val;y(last_d,:)];
    
    if active_stoch==1
    x_arr=double(int16((y_val(end,1:5).*[1 1/2 1/4 1/8 1/16])/sum(y_val(end,1:5))*Xi));
    if sum(x_arr.*[1 2 4 8 16])<Xi 
        x_arr(1)=x_arr(1)+1;
    end
    X=x_arr(1);X_2=x_arr(2);X_4=x_arr(3);X_8=x_arr(4);X_16=x_arr(5);
    end
    f_sedtotal=[f_sedtotal; (sum(sum(nssij_sed.*arr_sizes/arr_sizes(1)))+sum(n_sed.*arr_sizes/arr_sizes(1)))]; 
    f_total_np_system=[f_total_np_system; f_sedtotal(end)+f_hete(end)+f_single(end)+f_homo(end)]; 
   
    
C1=fg(1)*Z1;
C2=fg(2)*Z2;
C3=fg(3)*Z3;
C4=fg(4)*Z4;
C5=fg(2)*Z_arr(5);
C6=fg(1)*Z_arr(6);
C7=fg(1)*Z_arr(7);
C8=fg(1)*Z_arr(8);
C_arr=[C1;C2;C3;C4;C5;C6;C7;C8];
m_arr=C_arr.*V_arr;
Results_level3=[Results_level3; Vnp/Vwater m_arr' C_arr' fg(5)];    
    tt=tt+dt;
    %adj
end
avg_Knp=mean(Knp_process);
ff_Knp=[ff_Knp; avg_Knp];
ff=[ff (ones(numel(f_single),1)*inital_c-(f_hete+f_single+f_homo))];
ff_het=[ff_het f_hete];
ff_hom=[ff_hom f_homo];
ff_single=[ff_single f_single];
ff_water=[ff_water Results_level3(:,10)];
ff_air=[ff_air Results_level3(:,11)];
ff_sediment=[ff_sediment Results_level3(:,13)];
ff_sus_sediment=[ff_sus_sediment Results_level3(:,15)];
ff_fish=[ff_fish Results_level3(:,16)];
ff_np=[ff_np Results_level3(:,17)];
ff_sed=[ff_sed Results_level3(:,13)];
ff_fish_np=[ff_fish_np Results_level3(:,18).*Results_level3(:,17)];
ff_np_in_fish=[ff_np_in_fish f_fish_np];
ff_sedtotal=[ff_sedtotal f_sedtotal];
ff_total_np_system=[ff_total_np_system f_total_np_system];
ff_avg_size=[ff_avg_size f_avg_size];

%arr_K(ccc+1,:,:)=K;
%arr_Kss(ccc+1,:,:)=Kss;
ccc=ccc+1;
[plastic_num sal ccc]
end
   % solvedndt(K, n,alpha,dt,n_prev)
    %euler
    %v_dndt=dndt(K,n,alpha);
    %n=v_dndt'*dt+n_prev;  
% if n(numel(n))>1e-10
%    n=[n;0];
%    arr_sizes=[arr_sizes; size_p*(numel(n))^(1/3)];
%    K=dndt(K, n,alpha,arr_sizes);
% end

    %if mod(ii,(60/dt))==0
     %  hydro_size=[hydro_size; t (sum(arr_sizes.*n)*1e9)/sum(n)];
       
    %end
    
    %t=t+dt
    %ii=ii+1;
    
%f=[f_single(end) f_hete(end) f_homo(end) f_sedtotal(end) f_total_np_system(end) f_avg_size(end) alpha_ss]; 
    
%end
%plot_agg3;