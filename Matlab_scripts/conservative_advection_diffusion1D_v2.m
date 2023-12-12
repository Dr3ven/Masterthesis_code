clear, figure(1), clf
% Physics
Lx      = 1;
rho0    = 1;
Vx0     = -0.3*0.03;
D       = 2e-4;
% numerics
nx      = 200;
dx      = Lx/(nx);
nout    = 10;
dt      = min(dx.^2/D/2.1,dx/abs(Vx0)/1.1);
% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx-dx)/2:dx:(Lx-dx)/2;
rho     = rho0*exp(-300*(Xv+0).^2);
rho     = rho0*exp(-300*(Xc-0.3).^2);
rho2    = rho;
Vx      = Vx0*ones(1,nx+1);
LW = 0;
Vx      = linspace(Vx0,Vx0/2,length(Vx));
plot(Xc, rho)
% Action
for it = 1:5e4
    Vxc          = 0.5*(Vx(1:end-1)+Vx(2:end));
    rhoc         = 0.5*(rho(1:end-1)+rho(2:end));
    qx           = -D*diff(rho)/dx;
%     Fx           = rho.*av(Vx) + [av(qx(1:2)) av(qx) av(qx(end-1:end))];
    Fx           = rho.*av(Vx) + [qx(1) av(qx) qx(end)];
    Fxc          = av(Fx) + (- sign(Vx(2:end-1)) + LW.*(sign(Vx(2:end-1))-Vx(2:end-1)*dt/dx) ).*diff(Fx)/2 ;
    rho(2:end-1) = rho(2:end-1) - diff(Fxc)/dx*dt; 
    rho(end)     = rho(end) - (Vx(end)>0)*Vx(end)*(rho(end)-rho(end-1))/dx*dt;
    rho(1)       = rho(1  ) - (Vx(1  )<0)*Vx(1  )*(rho(2  )-rho(1    ))/dx*dt;
    qx2           = -D*diff(rho2)/dx;
    rho2         = rho2 - rho2.*diff(Vx)/dx*dt;
    rho2(2:end-1)= rho2(2:end-1) - diff(qx2)/dx*dt;
    rho2(2:end)  = rho2(2:end)   - (Vx0>0)*Vx(2:end-1).*diff(rho2)/dx*dt;
    rho2(1:end-1)= rho2(1:end-1) - (Vx0<0)*Vx(2:end-1).*diff(rho2)/dx*dt;
    if mod(it-1, nout) == 0
        plot(Xc, rho,'-',Xc, rho2,'--'), axis([-Lx/2 Lx/2 0 rho0]), drawnow
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


