function f=find_sedimentation(i,j,ai,aj,rho_p,rho_ps,rhow,mu)
f=2*9.81*((ai/2)^3+(aj/2)^3)^(2/3)*((rho_ps*(4/3)*pi*aj^3+rho_p*(4/3)*pi*ai^3)/(4/3*pi*aj^3+4/3*pi*ai^3)-rhow)/(9*mu);
if f<0 f=0; end
%((ai/a1)^df*4/3*pi*a1^3*rho_p+rho_w*(4/3*pi*ai^3-(ai/a1)^df*4/3*pi*a1^3))/(4/3*pi*ai^3);
end