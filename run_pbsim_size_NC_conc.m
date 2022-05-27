for plastic_num=1:1

for NC_count=41:50
%plastic_num=1;
fm_sal=zeros(4,51,2); w = warning ('off','all');
fug_sal=zeros(7,51,2); 
Knp_sal=[];
fm_size=zeros(51,2);

for sal=1:2
%plastic_num=1;
T=298.15;
Psed=2400;
NCsize_a=500;
NCconc_a=NC_count-1;
Shear_a=0;
density_array=[1050 2200 932 1450 858 1350]; 
H_array=[5.91E-21 3.6E-21 8E-21 7.56E-21 6e-21 5e-21];
%D_arr=[3.55e-14 3.98e-14 3.98e-14 ];
%PS       PTFE       PE      PVC       Polyprop  polyester
%1050     2200       932     1450      858       1350 
densityNP_a=density_array(plastic_num); %
densityNC_a=2650; %
densityNOM_a=1.33;
%NPsize_a=100.0;
active_stoch=0;
resus_rate_a=0.1*(res_count-1);
% n_dom=1;
% Xi=50;
% X=Xi; %lower is faster X_2=0;
% X_4=0;X_8=0;X_16=0;X_32=0;X_64=0;Xdom=X*n_dom; Xdomi=Xdom; nn_domfrac=0;
% alpha_np_np=alpha; alpha_dom_dom=0; alpha_dom_np=1;
% Rdom=40.0e-9/2; C=inital_c/(1.82e15); %mg/L

fractal_d_a=2.5;
%PS       PTFE     PMMA      PE      PVA      PVC       Polyprop  polyester
%5.91E-21 3.6E-21  3.44E-21  8E-21   1.1E-20  7.56E-21  6e-21     5e-21
Hamaker_A_a=H_array(plastic_num); % PS
Hamaker_Ass_a=(sqrt(1.5e-19))^2; % Kao
depth_a=85;
Vwater=1.63e12;
depth=depth_a;
Rmax_a=0.04; 
OM_size_a=0e-9*(active_stoch==0);
OM_packing_a=0.01; 
Temperature_a=298.15;
concentrationNP_a=1000; %ug/L Plastic
%sal=500;
salinity_a=10; % vary this to show effect
if sal>1
   %salinity_a=(sal-1)*100; 
   salinity_a=600;
end
NPsize_a=100; % vary this to show effect 
%salinity_a=600; 

%  [8.30598e-7 6.6661e-5 2.7998e-5 7.3194e-7 1.27856e-7 1.3332e-8]
chem_num=4;
chem_properties=[228.29 252.31 228.29 252.31 252.31 276.33;...
    254 163 158 175 217 278;...
    1.2e-7  6.67e-5  2.8e-5  5.74E-5  1.3e-8  2.55e-5;...
    2e-3    1.5e-3   9.4e-3  3.8e-3   8e-4   0.079;...
    0.5299297, 0.06657053, 1.2159, 0.04630553, 0.0591738, 0.03353858;...
    5.8 5.8 5.9 6.0 6.2 6.9;...
    0.41*10^5.8 0.41*10^5.8 0.41*10^5.9 4.04e5 0.41*10^6.2 0.41*10^6.9;...
    7.62 7.65 7.3 7.95 7.61 8.04;...
    151.75     44.97  227.88   10.55     34.15       912.8;...
    220        52     170      25        43          45;...
    0.025 0.049  0.017    0.025     0.049       0.05];
    
Molar_mass=chem_properties(1,chem_num);%252.13; %228.29 252.31 228.29 252.31 252.31 276.33
%temp=9.48;  
Melting_pt=chem_properties(2,chem_num); % 254 163 158 175 217 278
VP=chem_properties(3,chem_num);     %1.2e-7  6.67e-5  2.8e-5  5.74E-5,7.3e-7  1.3e-8  2.55e-5
CW=chem_properties(4,chem_num);   %2e-3    1.5e-3   9.4e-3  3.8e-3   8e-4   0.079
       
H=chem_properties(5,chem_num);  %0.5299297, 0.06657053, 1.2159, 0.04630553, 0.0591738, 0.03353858
             %CHR BBF BAA BAP BKF BGP
logKow=chem_properties(6,chem_num); %5.8 5.8 5.9 6.0 6.2 6.9
Koc=chem_properties(7,chem_num);
%Koc=chem_properties(8,chem_num);
F_ratio=exp(6.79*(1-(Melting_pt+273.15)/(T)));
SLVP=VP/F_ratio;
M=222.26;
           % CHR BBF  BAA BAP  BKF  BGP
Knp=chem_properties(8,chem_num);  %7.62 7.65 7.3 7.95 7.61 8.04
                    % chrysene  BBF    BAA      BAP       BKF         BGP
Cwi=chem_properties(9,chem_num);  %pg/L    151.75     44.97  227.88   10.55     34.15       912.8
Cwi_mol=Cwi*1e-9/Molar_mass; 
Cai=chem_properties(10,chem_num);     %pg/m3   220        52     170      25        43          45    
Cai_mol=Cai*1e-12/Molar_mass;
Csedi=chem_properties(11,chem_num); %ug/g              0.049  0.017    0.025     0.049       0.05?
Csedi_mol=Csedi*1e-6*1000*Psed/Molar_mass;
%


pb_sim4;%(NPsize_a,NCsize_a,NCconc_a,shear_a,densityNP_a,densityNC_a,...

numpts=102;
fm=[ff_hom(numpts,:)./ff_total_np_system(numpts,:); ff_het(numpts,:)./ff_total_np_system(numpts,:); ff_single(numpts,:)./ff_total_np_system(numpts,:); ff_sedtotal(numpts,:)./ff_total_np_system(numpts,:)];
%when running once: see variation with respect to time
%fug_res=[Results_level3(:,10) Results_level3(:,11) Results_level3(:,13) Results_level3(:,15) Results_level3(:,16) Results_level3(:,17)]; %water air sediments NC fish NP
%fug_res=fug_res';

%fug_res=[ff_water(end,:); ff_air(end,:); ff_sediment(end,:); ff_sus_sediment(end,:); ff_fish(end,:); ff_np(end,:); ff_mass_np(end,:)];
%densityNOM_a,salinity_a,fractal_d_a,Hamaker_A_a,Hamaker_Ass_a,depth_a,Rmax_a, OM_size_a, OM_packing_a, Temperature_a)
fm_sal(:,:,sal)=fm;
fm_size(:,sal)=ff_avg_size(numpts,:);
%fug_sal(:,:,sal)=fug_res;
Knp_sal=[Knp_sal ff_Knp];
end
 fm1=[fug_sal(:,:,1); fug_sal(:,:,2)];% fug_sal(:,:,3); fug_sal(:,:,4); fug_sal(:,:,5)];
%  fm1=[fm1; fug_sal(:,:,6); fug_sal(:,:,7); fug_sal(:,:,8); fug_sal(:,:,9); fug_sal(:,:,10)];
%  fm1=[fm1; fug_sal(:,:,11); fug_sal(:,:,12); fug_sal(:,:,13); fug_sal(:,:,14); fug_sal(:,:,15)];
%  fm1=[fm1; fug_sal(:,:,16); fug_sal(:,:,17); fug_sal(:,:,18); fug_sal(:,:,19); fug_sal(:,:,20)];
 string_filename=strcat('plastic_fug_nv_',num2str(plastic_num),'_',num2str(NC_count-1),'_mgL_NC_size_effect.csv');
 csvwrite(string_filename,fm1);
 fm1=[fm_sal(:,:,1); fm_sal(:,:,2)];%; fm_sal(:,:,3); fm_sal(:,:,4); fm_sal(:,:,5)];
%  fm1=[fm1; fm_sal(:,:,6); fm_sal(:,:,7); fm_sal(:,:,8); fm_sal(:,:,9); fm_sal(:,:,10)];
%  fm1=[fm1; fm_sal(:,:,11); fm_sal(:,:,12); fm_sal(:,:,13); fm_sal(:,:,14); fm_sal(:,:,15)];
%  fm1=[fm1; fm_sal(:,:,16); fm_sal(:,:,17); fm_sal(:,:,18); fm_sal(:,:,19); fm_sal(:,:,20)];
 string_filename=strcat('plastic_fate_nv_',num2str(plastic_num),'_',num2str(NC_count-1),'_mgL_NC_size_effect.csv');
 csvwrite(string_filename,fm1);
 string_filename=strcat('plastic_fate_nv_',num2str(plastic_num),'_',num2str(NC_count-1),'_mgL_NC_size_of_aggregate.csv');
 csvwrite(string_filename,fm_size);
 string_filename=strcat('plastic_fate_nv_',num2str(plastic_num),'_',num2str(NC_count-1),'_mgL_NC_Knp_rates.csv');
 csvwrite(string_filename,Knp_sal);
end
end