clear, figure(1), clf
% Physics
Lx      = 1;
rho0    = 1;
Vx0     = -1e-2;
D       = 2e-5;
% numerics
nx      = 300;
dx      = Lx/(nx);
nout    = 1;
dt      = min(dx.^2/D/2.5,dx/abs(Vx0)/2.5);
% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx+dx)/2:dx:(Lx+dx)/2;
rho     = rho0*exp(-300*(Xc+0).^2);
Vx      = Vx0*ones(1,nx+1);
%Vx      = linspace(Vx0/2,Vx0,length(Xc));
plot(Xc, rho)
% Action
for it = 1:5e2
    rho_old = rho;
    err = 1;
    iter = 0;
    while err > 1e-10
        qx           = -D*diff(rho)/dx;
        Frhox        = (Vx>0).*Vx.*rho(1:end-1) + (Vx<0).*Vx.*rho(2:end);
%         Frhox        = Vx.*av(rho);
        drhodt       = (rho-rho_old)/dt;
        rhoRes       = - drhodt(2:end-1) - diff(Frhox)/dx; 
        rho(2:end-1) = rho(2:end-1) + rhoRes*dt;
        rho(end)     = rho(end) - (Vx(end)>0)*Vx(end)*(rho(end)-rho(end-1))/dx*dt;
        rho(1)       = rho(1  ) - (Vx(1  )<0)*Vx(1  )*(rho(2  )-rho(1    ))/dx*dt;
%         rho(1) = -rho(2);
%         rho(end)= -rho(end-1);
        err = max(abs(rhoRes(:)));
    end
    if mod(it-1, nout) == 0
        plot(Xc, rho,'-'), axis([-Lx/2 Lx/2 0 rho0]), drawnow
    end
end
function A = av(B)
    A = 0.5*(B(2:end)+B(1:end-1));
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


