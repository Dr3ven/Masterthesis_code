clear, figure(1), clf
% Physics
Lx      = 1;
rho0    = 1;
Vx0     = -0.3;
beta    = 1;
tt      = 0.6/abs(Vx0);
% numerics
nx      = 15e2;                               % number of cells
dx      = Lx/(nx);
nout    = 10;

nt      = 2e6;
LW      = 0;                                 % set to 1 to choose Lax-Wendroff (spurious oscillations) and to 0 to choose upwind (diffusive), or to any intermediate value
% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx-dx)/2:dx:(Lx-dx)/2;
x2      = Xc;
rho     = rho0*exp(-300*(Xc-0.3).^2)+rho0;
P       = rho; 
% rho     = rho0*heaviside(Xc-0.3)+1;
% for i = 1:20
% rho(2:end-1) = rho(2:end-1) + 1/4.1*diff(rho,2,2); 
% end
rho2    = rho;
%rho     = rho0*ones(size(Xg));
Vx      = Vx0*ones(size(Xv));
rhoVx   = rho.*av(Vx);
%Vx      = Vx0*linspace(1,0.5,nx);
Vxini   = Vx;
plot(Xc, rho, Xc, rhoVx)
% Action
for it = 1:nt
    dt             = min(dx/max(1./beta./abs(Vx))/1.1,dx/max(abs(Vx))/2.1);
    Frho           = rhoVx;
    Fcrho          = av(Frho) + (- sign(Vx(2:end-1)) + LW.*(sign(Vx(2:end-1))-Vx(2:end-1)*dt/dx) ).*diff(Frho)/2;
    rho(2:end-1)   = rho(2:end-1) - diff(Fcrho)/dx*dt;
    rho(end)       = rho(end) - (Vx(1)>0)*(Frho(end)-Frho(end-1))/dx*dt;
    rho(1)         = rho(1)   - (Vx(end)<0)*(Frho(2)-Frho(1))/dx*dt;
    dPdt           = -1/beta*diff(Vx)/dx;
    P              = P + dPdt*dt;
    FVx            = rhoVx.*av(Vx);
    FcVx           = av(FVx) + (- sign(Vx(2:end-1)) + LW.*(sign(Vx(2:end-1))-Vx(2:end-1)*dt/dx) ).*diff(FVx)/2;
    rhoVx(2:end-1) = rhoVx(2:end-1) - diff(FcVx)/dx*dt - diff(av(rho))/dx*dt;
%     rhoVx(end)     = rhoVx(end) - (Vx0>0)*(FVx(end)-FVx(end-1))/dx*dt;
    rhoVx(end)     = Vx0 + -sin(30*it*dt)*Vx0/10;
    rhoVx(1)       = rhoVx(1)   - (Vx0<0)*(FVx(2)-FVx(1))/dx*dt;
    Vx(2:end-1)    = av(rhoVx./rho);
    Vx(1)          = 2*rhoVx(1  )/rho(1  )-Vx(2    );
    Vx(end)        = 2*rhoVx(end)/rho(end)-Vx(end-1);
    %FcVx         = 0.5.*(rhoVx(2:end)+rhoVx(1:end-1)) + (- sign(Vx(2:end-1)) + LW.*(sign(Vx(2:end-1))-Vx(2:end-1)*dt/dx) ).*diff(rhoVx)/2;
    

    
    
    if mod(it-1, nout) == 0
%         plot(Xc, rho,'-k',Xc, rho2,'--r',Xc, rho0*heaviside(Xc-Vx0*dt*it)+1,'--b'), axis([-Lx/2 Lx/2 0 rho0+1.5 ]), drawnow
        plot(Xc, rho, Xc, rhoVx,Xc, P), drawnow
        %plot(phirho), drawnow
    end
end

[sum(rho*dx) sum(rho2*dx) sum((rho0*exp(-3000*(Xc-Vx0*dt*it-0.3).^2)+1)*dx)]-sum((rho0*exp(-3000*(Xc-Vx0*dt*it-0.3).^2)+1)*dx)

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


