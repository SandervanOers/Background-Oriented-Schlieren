clc;clear;close all;

t = 0:1:80;
t=deg2rad(t);
n=1.333;
Dx1 = tan(t) - sin(t)./(n*sqrt(1 - sin(t).^2/n^2));
n=1.4;
Dx2 = tan(t) - sin(t)./(n*sqrt(1 - sin(t).^2/n^2));

figure;
hold all
box on
plot(rad2deg(t),Dx1)
plot(rad2deg(t),Dx2)
xlabel('$\theta_x [^{\circ}]$', 'interpreter', 'latex')
ylabel('$\frac{\Delta x}{L_t}$', 'interpreter', 'latex')
legend('n = 1.333', 'n = 1.4')
set(gca,'FontSize',40)
hold off

%%
Dx = Dx1;
dndt1 = -(3 + (Dx.^2 - 2).*cos(t).^2 - 3*Dx.*cos(t).*sin(t)).*cos(t).*Dx./(sqrt((Dx.^2.*cos(t).^2 - 2.*Dx.*cos(t).*sin(t) + 1)./cos(t).^2).*(-2.*Dx.*cos(t).*sin(t) + 1 + (Dx.^2 - 1).*cos(t).^2));
Dx = Dx2;
dndt2 = -(3 + (Dx.^2 - 2).*cos(t).^2 - 3*Dx.*cos(t).*sin(t)).*cos(t).*Dx./(sqrt((Dx.^2.*cos(t).^2 - 2.*Dx.*cos(t).*sin(t) + 1)./cos(t).^2).*(-2.*Dx.*cos(t).*sin(t) + 1 + (Dx.^2 - 1).*cos(t).^2));

figure;
hold all
box on
plot(rad2deg(t),dndt1, 'k')
plot(rad2deg(t),dndt2, 'k--')
xlabel('$\theta_x [^{\circ}]$', 'interpreter', 'latex')
ylabel('$\frac{\partial n}{\partial \theta_x}$', 'interpreter', 'latex')
legend('n = 1.333', 'n = 1.4')
set(gca,'FontSize',40)
hold off

%%
Dx = Dx1;
dnddx1=sin(t).*cos(t).^2./(sqrt((Dx.^2.*cos(t).^2 - 2*Dx.*cos(t).*sin(t) + 1)./cos(t).^2).*(-2*Dx.*cos(t).*sin(t) + 1 + (Dx.^2 - 1).*cos(t).^2));
Dx = Dx2;
dnddx2=sin(t).*cos(t).^2./(sqrt((Dx.^2.*cos(t).^2 - 2*Dx.*cos(t).*sin(t) + 1)./cos(t).^2).*(-2*Dx.*cos(t).*sin(t) + 1 + (Dx.^2 - 1).*cos(t).^2));

figure;
hold all
box on
plot(rad2deg(t),dnddx1, 'k')
plot(rad2deg(t),dnddx2, 'k--')
xlabel('$\theta_x [^{\circ}]$', 'interpreter', 'latex')
ylabel('$\frac{\partial n}{\partial (\Delta x / L_t)}$', 'interpreter', 'latex')
legend('n = 1.333', 'n = 1.4')
set(gca,'FontSize',40)
hold off