clear, figure(1), clf
% Physics
Lx      = 1;
rho0    = 1;
Vx0     = 1e-2;
% numerics
nx      = 1000;
dx      = Lx/(nx);
nout    = 3;
dt      = 1e-1;
% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx-dx)/2:dx:(Lx-dx)/2;
rho     = rho0*exp(-300*(Xc+0).^2);
intrho0 = sum(rho*dx);
Vx      = Vx0*ones(1,nx+1);
Vx      = linspace(Vx0/2,Vx0,length(Xv));
rho2    = rho;
plot(Xc, rho)
% Action
for it = 1:5e4
    Vxc         = 0.5*(Vx(1:end-1)+Vx(2:end));
    rhoc        = 0.5*(rho(1:end-1)+rho(2:end));
    rho2c       = 0.5*(rho2(1:end-1)+rho2(2:end));
    drhodt      = -diff(Vx(2:end-1).*rhoc)/dx;
    rho(2:end-1)  = rho(2:end-1) + drhodt*dt;
    drhodt2     = -Vx(2:end-1).*diff(rho2)/dx - rho2c.*diff(Vxc)/dx;
    rho2(2:end) = rho2(2:end) + drhodt2*dt;
    subplot(211)
    plot(Xc, rho,'-', Xc, rho2,'x'), axis([-Lx/2 Lx/2 0 rho0])
    subplot(212)
    plot(it, sum(rho*dx)-intrho0, 'or',it, sum(rho2*dx)-intrho0, 'xb'), hold on, drawnow
end


function [A] = advection2d(A,Vx,Vy,dx,dy,dt)
%Function advection 2d is a function that is used for upwind
%advection in 2d problems
dtadv = dt;
A(1:end-1,:) = A(1:end-1,:) - (Vx(2:end-1,:)<0).*Vx(2:end-1,:).*diff(A,1,1)/dx*dtadv;
A(2:end  ,:) = A(2:end  ,:) - (Vx(2:end-1,:)>0).*Vx(2:end-1,:).*diff(A,1,1)/dx*dtadv;        
A(:,1:end-1) = A(:,1:end-1) - (Vy(:,2:end-1)<0).*Vy(:,2:end-1).*diff(A,1,2)/dy*dtadv;
A(:,2:end  ) = A(:,2:end  ) - (Vy(:,2:end-1)>0).*Vy(:,2:end-1).*diff(A,1,2)/dy*dtadv;
return
end


