%% validate_output.m — Run in MATLAB, paste output back
%  Produces numerical values for all key integrals for comparison with Julia.

clear; clc;

%% Parameters
U=26.7046; g=981.0; T=72.0; rho_l=1.0; rho_u=0.001;
l_c=U^2/g; t_c=U/g; F_c=rho_l*U^2*l_c;
alpha=T/(rho_l*U^2*l_c); rho_r=rho_u/rho_l;
beta=(1-rho_r)/(1+rho_r); gamma_rho=1/(1+rho_r);
sqrt_beta=sqrt(beta);
disc=(1+rho_r)^2-4*alpha*(1-rho_r);
k_l=((1+rho_r)+sqrt(disc))/(2*alpha);
k_s=((1+rho_r)-sqrt(disc))/(2*alpha);
F0=0.01*T/F_c;

fprintf('=== CG PARAMETERS ===\n');
fprintf('alpha   = %.15e\n', alpha);
fprintf('k_s     = %.15e\n', k_s);
fprintf('k_l     = %.15e\n', k_l);
fprintf('F0      = %.15e\n', F0);

%% chi(2)
chi=@(k) sqrt(beta.*k+gamma_rho.*alpha.*k.^3);
fprintf('\nchi(2)  = %.15e\n', chi(2.0));

%% Combined integrand at k=2, x=3, t=110
x=3.0; t=110.0; epsilon_pv=1e-6;
total_integrand=@(k) ...
  2.0.*cos(k.*x)./(alpha.*(k-k_l).*(k-k_s)) ...
  -(1+rho_r).*(k+chi(k))./(1-rho_r+alpha.*k.^2).*cos(k.*(t-x)-t.*chi(k))./(alpha.*(k-k_l).*(k-k_s)) ...
  -(1+rho_r).*(k-chi(k))./(1-rho_r+alpha.*k.^2).*cos(k.*(t-x)+t.*chi(k))./(alpha.*(k-k_l).*(k-k_s));

fprintf('\nI(2; x=3, t=110) = %.15e\n', total_integrand(2.0));

%% Partial integrals
I1=integral(total_integrand, 0, k_s-epsilon_pv, 'AbsTol',1e-10,'RelTol',1e-8);
I2=integral(total_integrand, k_s+epsilon_pv, k_l-epsilon_pv, 'AbsTol',1e-10,'RelTol',1e-8);
I3=integral(total_integrand, k_l+epsilon_pv, Inf, 'AbsTol',1e-10,'RelTol',1e-8);
eta_ivp=-F0/(2*pi)*(I1+I2+I3);
fprintf('\n=== CG PARTIAL INTEGRALS (x=3, t=110) ===\n');
fprintf('I1      = %.15e\n', I1);
fprintf('I2      = %.15e\n', I2);
fprintf('I3      = %.15e\n', I3);
fprintf('eta_IVP = %.15e\n', eta_ivp);

%% Steady
G_integrand=@(k) cos(k.*x)./(k+k_s)-cos(k.*x)./(k+k_l);
G_x=1/(k_l-k_s)*integral(G_integrand, 0, Inf, 'AbsTol',1e-10,'RelTol',1e-8);
eta_steady=-2*F0/(alpha*(k_l-k_s))*sin(k_s*x)+F0/(pi*alpha)*G_x;
fprintf('G(3)        = %.15e\n', G_x);
fprintf('eta_steady  = %.15e\n', eta_steady);

%% ========== PURE GRAVITY ==========
fprintf('\n=== PG PARAMETERS ===\n');
fprintf('beta       = %.15e\n', beta);
fprintf('sqrt_beta  = %.15e\n', sqrt_beta);

%% T terms at x=-2, t_dim=1.0s
t_pg=1.0/t_c; x_pg=-2.0; a_pg=t_pg-x_pg;
fprintf('\n=== PG T-TERMS (x=-2, t_dim=1s, t=%.6f, a=%.6f) ===\n', t_pg, a_pg);

T0a=-pi*sin(beta*abs(x_pg));
T0b=integral(@(y) exp(-y.*abs(x_pg)).*y./(beta^2+y.^2), 0, Inf, 'AbsTol',1e-10,'RelTol',1e-8);
T0=(T0a+T0b)/(pi*(1+rho_r));
fprintf('T0  = %.15e\n', T0);

T1=-sin(beta*x_pg)/(1+rho_r);
fprintf('T1- = %.15e\n', T1);

T2=-4/(pi*(1+rho_r)*beta)*integral(@(v) ...
  v.^2.*exp(-2*v.^2*a_pg+v*t_pg*sqrt_beta).* ...
  (sqrt_beta*cos(v*t_pg*sqrt_beta)+(2*v-sqrt_beta).*sin(v*t_pg*sqrt_beta)) ...
  ./(beta+(2*v-sqrt_beta).^2), 0, 100, 'AbsTol',1e-10,'RelTol',1e-8);
fprintf('T2- = %.15e\n', T2);

b=0.5*t_pg*sqrt_beta; X=b*sqrt(2/(pi*a_pg));
T3=1/(pi*(1+rho_r)*sqrt_beta)*(1+t_pg/(2*a_pg))*sqrt(pi/(2*a_pg))* ...
  (cos(b^2/a_pg)*(0.5-fresnelc(X))+sin(b^2/a_pg)*(0.5-fresnels(X)));
fprintf('T3- = %.15e\n', T3);

T4=-1/(pi*(1+rho_r))*integral(@(v) ...
  cos(v.^2*a_pg+v*t_pg*sqrt_beta)./(v+sqrt_beta), 0, 100, 'AbsTol',1e-10,'RelTol',1e-8);
fprintf('T4- = %.15e\n', T4);

eta_ana=F0*(T0+T1+T2+T3+T4);
fprintf('eta_analytical = %.15e\n', eta_ana);

%% CPV
combined=@(k) cos(k.*x_pg)./(pi*(1+rho_r).*(k-beta)) ...
  -k.*cos(k.*(t_pg-x_pg)-t_pg.*sqrt(beta.*k))./(2*pi*(1-rho_r).*(k-sqrt(beta.*k))) ...
  -k.*cos(k.*(t_pg-x_pg)+t_pg.*sqrt(beta.*k))./(2*pi*(1-rho_r).*(k+sqrt(beta.*k)));
I_below=integral(combined, 0, beta-1e-6, 'AbsTol',1e-10,'RelTol',1e-8);
I_above=integral(combined, beta+1e-6, 100, 'AbsTol',1e-10,'RelTol',1e-8);
eta_cpv=F0*(I_below+I_above);
fprintf('eta_CPV        = %.15e\n', eta_cpv);
fprintf('|diff|          = %.4e\n', abs(eta_ana-eta_cpv));

fprintf('\nDone.\n');
