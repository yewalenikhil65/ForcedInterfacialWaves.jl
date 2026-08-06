%% generate_comparison_data.m
%  Generates profile data at two time indices and saves as CSV.
%  Julia will plot both on the same axes.
%  Run in MATLAB, then copy the .csv files to docs/src/assets/

clear; clc; close all;

U=26.7046; g=981.0; T=72.0; rho_l=1.0; rho_u=0.001;
l_c=U^2/g; t_c=U/g; F_c=rho_l*U^2*l_c;
alpha=T/(rho_l*U^2*l_c); rho_r=rho_u/rho_l;
beta=(1-rho_r)/(1+rho_r); gamma_rho=1/(1+rho_r);
sqrt_beta=sqrt(beta);
disc=(1+rho_r)^2-4*alpha*(1-rho_r);
k_l=((1+rho_r)+sqrt(disc))/(2*alpha);
k_s=((1+rho_r)-sqrt(disc))/(2*alpha);
F0=0.01*T/F_c; epsilon_pv=1e-6;

chi=@(k) sqrt(beta.*k+gamma_rho.*alpha.*k.^3);
total_integrand=@(k,x,t) ...
  2.0.*cos(k.*x)./(alpha.*(k-k_l).*(k-k_s)) ...
  -(1+rho_r).*(k+chi(k))./(1-rho_r+alpha.*k.^2).*cos(k.*(t-x)-t.*chi(k))./(alpha.*(k-k_l).*(k-k_s)) ...
  -(1+rho_r).*(k-chi(k))./(1-rho_r+alpha.*k.^2).*cos(k.*(t-x)+t.*chi(k))./(alpha.*(k-k_l).*(k-k_s));

Nx=201;
x_grid=linspace(-5, 10, Nx);

%% CG Figure 10 at time_index=300
t_dim=3.0; t=t_dim/t_c;
eta_ivp=zeros(1,Nx);
for ix=1:Nx
  f=@(k) total_integrand(k, x_grid(ix), t);
  I1=integral(f, 0, k_s-epsilon_pv, 'AbsTol',1e-10,'RelTol',1e-8);
  I2=integral(f, k_s+epsilon_pv, k_l-epsilon_pv, 'AbsTol',1e-10,'RelTol',1e-8);
  I3=integral(f, k_l+epsilon_pv, Inf, 'AbsTol',1e-10,'RelTol',1e-8);
  eta_ivp(ix)=-F0/(2*pi)*(I1+I2+I3);
end
csvwrite('matlab_cg_ivp_t300.csv', [x_grid', eta_ivp']);
fprintf('Saved matlab_cg_ivp_t300.csv (%d points)\n', Nx);

%% PG Figure 6 at t_dim=1.0s (numerical CPV only, for overlay)
t_pg=1.0/t_c;
L=4*2*pi/beta;
x_pg=linspace(-L/2, L/2, Nx);
eta_cpv_pg=nan(1,Nx);
combined=@(k,xx) cos(k.*xx)./(pi*(1+rho_r).*(k-beta)) ...
  -k.*cos(k.*(t_pg-xx)-t_pg.*sqrt(beta.*k))./(2*pi*(1-rho_r).*(k-sqrt(beta.*k))) ...
  -k.*cos(k.*(t_pg-xx)+t_pg.*sqrt(beta.*k))./(2*pi*(1-rho_r).*(k+sqrt(beta.*k)));
for ix=1:Nx
  a=t_pg-x_pg(ix);
  if abs(a)>1.0
    f=@(k) combined(k, x_pg(ix));
    I_lo=integral(f, 0, beta-1e-6, 'AbsTol',1e-10,'RelTol',1e-8);
    I_hi=integral(f, beta+1e-6, 100, 'AbsTol',1e-10,'RelTol',1e-8);
    eta_cpv_pg(ix)=F0*(I_lo+I_hi);
  end
end
csvwrite('matlab_pg_cpv_t1s.csv', [x_pg', eta_cpv_pg']);
fprintf('Saved matlab_pg_cpv_t1s.csv (%d points)\n', Nx);

fprintf('Done. Copy CSV files to docs/src/assets/\n');
