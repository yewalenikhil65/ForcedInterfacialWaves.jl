%% Fig. 8(ii) — Julia vs MATLAB capillary-gravity IVP overlay
clear; close all; clc;
J = readmatrix('julia_fig8_profiles.csv');
xj = J(:,1); etaj = J(:,2); etatrj = J(:,3);
U=26.7046; g=981; T=72; rho_l=1; rho_u=.001;
lc=U^2/g; a=T/(rho_l*U^2*lc); r=rho_u/rho_l;
b=(1-r)/(1+r); gr=1/(1+r); F=.01*T/(rho_l*U^2*lc);
D=(1+r)^2-4*a*(1-r); kl=((1+r)+sqrt(D))/(2*a); ks=((1+r)-sqrt(D))/(2*a);
eps=1e-6; t=367.35; x=linspace(-15,15,2001); x(abs(x)<1e-12)=[];
assert(max(abs(x(:)-xj(:)))<1e-12);
chi=@(k)sqrt(b*k+gr*a*k.^3);
C=@(k)2*cos(k*x)./(a*(k-kl).*(k-ks)) ...
 -(1+r)*(k+chi(k)).*cos(k*(t-x)-t*chi(k))./((1-r+a*k^2)*a.*(k-kl).*(k-ks)) ...
 -(1+r)*(k-chi(k)).*cos(k*(t-x)+t*chi(k))./((1-r+a*k^2)*a.*(k-kl).*(k-ks));
I1=integral(C,0,ks-eps,'ArrayValued',true,'AbsTol',1e-10,'RelTol',1e-8);
I2=integral(C,ks+eps,kl-eps,'ArrayValued',true,'AbsTol',1e-10,'RelTol',1e-8);
I3=integral(C,kl+eps,Inf,'ArrayValued',true,'AbsTol',1e-10,'RelTol',1e-8);
eta=-F/(2*pi)*(I1+I2+I3);
G=@(k)(cos(k*x)./(k+ks)-cos(k*x)./(k+kl))/(kl-ks);
g=integral(G,0,Inf,'ArrayValued',true,'AbsTol',1e-10,'RelTol',1e-8);
etas=F/(a*(kl-ks))*(-sin(ks*abs(x))+sin(kl*abs(x)))+F*g/(pi*a);
etatr=eta-etas;
fprintf('max |eta_J-M| = %.6e\n',max(abs(etaj(:)-eta(:))));
fprintf('max |eta_tr_J-M| = %.6e\n',max(abs(etatrj(:)-etatr(:))));
idx=1:25:numel(x); idxTr=1:12:numel(x); figure('Color','w','Position',[100 100 760 360]); hold on; box on;
h1=plot(xj,etaj*1e3,'b--','LineWidth',2.5); h2=plot(xj,etatrj*1e3,'m:','LineWidth',2.5);
h3=plot(x(idx),eta(idx)*1e3,'o','LineStyle','none','Color',[0 0.35 0.95], ...
        'MarkerSize',5,'MarkerFaceColor','none','LineWidth',1.0);
h4=plot(x(idxTr),etatr(idxTr)*1e3,'s','LineStyle','none','Color',[0.90 0.10 0.10], ...
        'MarkerSize',5,'MarkerFaceColor','none','LineWidth',1.0);
xlim([-10 10]); ylim([-4.8 8.2]); yticks([-4 0 4 8]); set(gca,'FontSize',18,'TickLabelInterpreter','latex');
xlabel('$x$','Interpreter','latex','FontSize',22); ylabel('$\eta \times 10^3$','Interpreter','latex','FontSize',22);
legend([h1 h2 h3 h4],{'$\eta$ (Julia)','$\eta_{tr}$ (Julia)','$\eta$ (MATLAB)','$\eta_{tr}$ (MATLAB)'},'Interpreter','latex','FontSize',16,'Location','eastoutside','Box','off');
exportgraphics(gcf,'../fig8_overlay.png','Resolution',300); savefig(gcf,'fig8_overlay.fig');
