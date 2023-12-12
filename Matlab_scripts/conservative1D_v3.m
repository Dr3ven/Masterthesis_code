clear, figure(1), clf
% Physics
Lx      = 1;
rho0    = 1;
Vx0     = 1e-2;
D       = 2e-5;
beta    = 1e-3;
mu      = 1e-2;
% numerics
nx      = 300;
dx      = Lx/(nx);
nout    = 25;
c       = 1/sqrt(rho0*beta);
dt      = min([dx.^2/D/2.1,dx/abs(Vx0)/2.1,dx/c/2.1])/2;

% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx+dx)/2:dx:(Lx+dx)/2;
rho     = 0.5*exp(-300*(Xc+0).^2)+rho0;
P       = 1/beta*log(rho(2:end-1)/rho0);
Vx      = 0*Vx0*ones(1,nx+1);
%Vx      = linspace(Vx0/2,Vx0,length(Xc));
plot(Xc, rho, Xc(2:end-1), P)
% Action
for it = 1:5e4

    qx           = -D*diff(rho)/dx;
%     Frhox        = (Vx>0).*Vx.*rho(1:end-1) + (Vx<0).*Vx.*rho(2:end);
    Frhox        = Vx.*av(rho);
    drhodt       = -diff(Frhox + qx)/dx;
    rho(2:end-1) = rho(2:end-1) + drhodt*dt;    

%     rho(end)     = rho(end) - (Vx(end)>0)*Vx(end)*(rho(end)-rho(end-1))/dx*dt;
%     rho(1)       = rho(1  ) - (Vx(1  )<0)*Vx(1  )*(rho(2  )-rho(1    ))/dx*dt;
    rho(1)       = rho(2);
    rho(end)     = rho(end-1);
    rhoc         = 0.5*(rho(1:end-1)+rho(2:end));
    Vxc          = 0.5*(Vx(1:end-1)+Vx(2:end));
    Mx           = rhoc.*Vx;
    P            = 1/beta*log(rho(2:end-1)/rho0);
    Sxx          = -P + mu.*diff(Vx/dx);
%     FMxx         = (Vxc>0).*Vxc.*Mx(1:end-1) + (Vxc<0).*Vxc.*Mx(2:end);
    FMxx         =  Vxc.*av(Mx);
    dMxdt        = -diff(FMxx-Sxx)/dx;

    Mx(2:end-1)  = Mx(2:end-1) + dMxdt*dt;

%     Mx(end)      = Mx(end) - (Vxc(end)>0)*Vxc(end)*(Mx(end)-Mx(end-1))/dx*dt;
%     Mx(1)        = Mx(1  ) - (Vxc(1  )<0)*Vxc(1  )*(Mx(2  )-Mx(1    ))/dx*dt;
    Vx           = Mx./rhoc;
    if mod(it-1, nout) == 0
        plot(Xc, rho,'-',Xv, Mx,'-'), axis([-Lx/2 Lx/2 -2*rho0 2*rho0]), drawnow
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


