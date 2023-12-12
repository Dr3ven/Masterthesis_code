clear, figure(1), clf
% Physics
Lx      = 1;
rho0    = 1;
Vx0     = -0.3*0.03;
beta    = 1e-3; 
D       = 2e-4;
% numerics
nx      = 200;
dx      = Lx/(nx);
nout    = 1;

% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx-dx)/2:dx:(Lx-dx)/2;
Xg      = -(Lx+dx)/2:dx:(Lx+dx)/2;
P       = 1e2*exp(-300*Xg.^2);
rho     = rho0*(exp(beta.*P));
% rho     = rho0*exp(-300*(Xv+0).^2);
% rho     = rho0*exp(-300*(Xg-0.3).^2);
% rho2    = rho;
Vx      = Vx0*ones(1,nx+1);
LW = 0;
Vx      = linspace(0*Vx0,0*Vx0/2,length(Vx));
plot(Xg, rho)
% Action
dt      = min(dx/abs(Vx0)/1.1,dx./sqrt(1./(rho0*beta))/1.1)/10;
for it = 1:5e4
    Vxc           = 0.5*(Vx(1:end-1)+Vx(2:end));
    rhoc          = 0.5*(rho(1:end-1)+rho(2:end));
    qx            = -D*diff(rho)/dx;
%     Fx           = rho.*av(Vx) + [av(qx(1:2)) av(qx) av(qx(end-1:end))];
    Fx            = rhoc.*Vx + 0*qx;
    Fxc           = av(Fx) + (- sign(av(Vx)) + LW.*(sign(av(Vx))-av(Vx)*dt/dx) ).*diff(Fx)/2 ;
    rhoc(2:end-1) = rhoc(2:end-1) - diff(Fxc)/dx*dt; 
    rhoc(end)     = rhoc(end) - (Vx(end)>0)*Vx(end)*(rhoc(end)-rhoc(end-1))/dx*dt;
    rhoc(1)       = rhoc(1  ) - (Vx(1  )<0)*Vx(1  )*(rhoc(2  )-rhoc(1    ))/dx*dt;
    rho(2:end-1)  = av(rhoc);
    rho(end)      = rho(end) - (Vx(end)>0)*Vx(end)*(rho(end)-rho(end-1))/dx*dt;
    rho(1)        = rho(1  ) - (Vx(1  )<0)*Vx(1  )*(rho(2  )-rho(1    ))/dx*dt;
    Mx            = rhoc.*Vx;
%     P             = 1/beta*log(rho/rho0);
    dPdt          = 1/beta*diff(Vx)/dx;
    P(2:end-1)    = P(2:end-1) + dPdt*dt;
    Fx            = rhoc.*Vx.*Vx + av(P);
    Fxc           = av(Fx) + (- sign(av(Vx)) + LW.*(sign(av(Vx))-av(Vx)*dt/dx) ).*diff(Fx)/2 ;
    Mx(2:end-1)   = Mx(2:end-1) - diff(Fxc)/dx*dt;
    Mx(end)       = Mx(end) - (Vx(end)>0)*Vx(end)*(Mx(end)-Mx(end-1))/dx*dt;
    Mx(1)         = Mx(1  ) - (Vx(1  )<0)*Vx(1  )*(Mx(2  )-Mx(1    ))/dx*dt;
    Vx            = Mx./rhoc;
    if mod(it-1, nout) == 0
        plot(Xv, av(P),'-',Xv, Mx,'--'), drawnow;%axis([-Lx/2 Lx/2 0 rho0]), drawnow
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


