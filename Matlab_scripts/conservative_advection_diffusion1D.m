clear, figure(1), clf
% Physics
Lx      = 1;
rho0    = 1;
Vx0     = -0.0003;
D       = 0*2e-4;
% numerics
nx      = 300;
dx      = Lx/(nx);
nout    = 10;
dt      = min(dx.^2/D/2.5,dx/abs(Vx0)/200.5);
% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx-dx)/2:dx:(Lx-dx)/2;
rho     = rho0*exp(-300*(Xv+0).^2);
rho     = rho0*exp(-300*(Xv-0.3).^2);
rho     = heaviside(Xv/Lx); 
rho2    = rho;
Vx      = Vx0*ones(1,nx);
%Vx      = linspace(Vx0/2,Vx0,length(Xc));
plot(Xv, rho)
% Action
for it = 1:5e2
    Vxc          = 0.5*(Vx(1:end-1)+Vx(2:end));
    rhoc         = 0.5*(rho(1:end-1)+rho(2:end));
    qx           = -D*diff(rho)/dx;
    drhodt       = -diff( (Vx>0).*Vx.*rho(1:end-1) + (Vx<0).*Vx.*rho(2:end) + qx)/dx;
    rho(2:end-1) = rho(2:end-1) + drhodt*dt;
    rho(end)     = rho(end) - (Vx(end)>0)*Vx(end)*(rho(end)-rho(end-1))/dx*dt;
    rho(1)       = rho(1  ) - (Vx(1  )<0)*Vx(1  )*(rho(2  )-rho(1    ))/dx*dt;
    rho2(2:end)  = rho2(2:end)   - (Vx0>0)*Vx0*diff(rho2)/dx*dt;
    rho2(1:end-1)= rho2(1:end-1) - (Vx0<0)*Vx0*diff(rho2)/dx*dt;
    if mod(it-1, nout) == 0
        plot(Xv, rho,'-',Xv, rho2,'--'), axis([-Lx/2 Lx/2 0 rho0]), drawnow
    end
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


