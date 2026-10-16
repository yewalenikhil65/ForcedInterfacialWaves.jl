%% Fig. 10(ii) — Julia vs MATLAB vs Basilisk overlay at t_dim = 0.25 s
clear; close all; clc;
J = readmatrix('julia_fig10_profiles.csv');
xj = J(:,1); etaj = J(:,2);
U=26.7046; g=981; T=72; rho_l=1; rho_u=.001;
lc=U^2/g; tc=U/g; a=T/(rho_l*U^2*lc); r=rho_u/rho_l;
b=(1-r)/(1+r); gr=1/(1+r); F=.01*T/(rho_l*U^2*lc);
D=(1+r)^2-4*a*(1-r); kl=((1+r)+sqrt(D))/(2*a); ks=((1+r)-sqrt(D))/(2*a);
eps=1e-6; tdim=0.25; t=tdim/tc;
x=linspace(-15,15,2001); x(abs(x)<1e-12)=[];
chi=@(k)sqrt(b*k+gr*a*k.^3);
C=@(k)2*cos(k*x)./(a*(k-kl).*(k-ks)) ...
 -(1+r)*(k+chi(k)).*cos(k*(t-x)-t*chi(k))./((1-r+a*k^2)*a.*(k-kl).*(k-ks)) ...
 -(1+r)*(k-chi(k)).*cos(k*(t-x)+t*chi(k))./((1-r+a*k^2)*a.*(k-kl).*(k-ks));
I1=integral(C,0,ks-eps,'ArrayValued',true,'AbsTol',1e-10,'RelTol',1e-8);
I2=integral(C,ks+eps,kl-eps,'ArrayValued',true,'AbsTol',1e-10,'RelTol',1e-8);
I3=integral(C,kl+eps,Inf,'ArrayValued',true,'AbsTol',1e-10,'RelTol',1e-8);
eta_m=-F/(2*pi)*(I1+I2+I3);
fprintf('max |eta_Julia - eta_MATLAB| = %.6e\n', max(abs(etaj(:)-eta_m(:))));

data=readmatrix(fullfile('..','..','..','..','notebooks','if_25.csv'));
data=sortrows(data,6); xb=data(:,6)/lc; yb=data(:,7)/lc;

idx=1:25:numel(x);
figure('Color','w','Position',[100 100 760 360]); hold on; box on;
h1=plot(xj,etaj*1e3,'b-','LineWidth',2.5);
h2=plot(xb,yb*1e3,'r:','LineWidth',2.5);
h3=plot(x(idx),eta_m(idx)*1e3,'o','LineStyle','none','Color',[0 0.35 0.95], ...
        'MarkerSize',5,'MarkerFaceColor','none','LineWidth',1.0);
xlim([-6 10]); ylim([-6 8.2]); yticks([-4 0 4 8]);
set(gca,'FontSize',18,'TickLabelInterpreter','latex');
xlabel('$x$','Interpreter','latex','FontSize',22);
ylabel('$\eta \times 10^3$','Interpreter','latex','FontSize',22);
legend([h1 h2 h3],{'$\eta$ (Julia)','Simulation (Basilisk)','$\eta$ (MATLAB)'}, ...
       'Interpreter','latex','FontSize',16,'Location','northeast','Box','off');
exportgraphics(gcf,'../fig10_overlay.png','Resolution',300);
savefig(gcf,'fig10_overlay.fig');
disp('Saved fig10_overlay.png and fig10_overlay.fig');
