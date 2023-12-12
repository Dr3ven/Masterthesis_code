clear, figure(1), clf
% Physics
Lx      = 1;
rho0    = 1;
Vx0     = -0.3;
D       = 2e-5;
beta    = 1;
eta     = 1;
tt      = 0.6/abs(Vx0);
% numerics
nx      = 5e2;                               % number of cells
dx      = Lx/(nx);
nout    = 10;
dt      = min(dx.^2/D/2.5,dx/abs(Vx0)/2.1);
nt      = tt/dt; 
% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx-dx)/2:dx:(Lx-dx)/2;
x2 = Xc;
rho     = rho0*exp(-3000*(Xc-0.3).^2)+1;

% rho     = rho0*heaviside(Xc-0.3)+1;
% for i = 1:20
% rho(2:end-1) = rho(2:end-1) + 1/4.1*diff(rho,2,2); 
% end
rho2    = rho;
%rho     = rho0*ones(size(Xg));
Vx      = Vx0*ones(size(Xv));
Vxini   = Vx;
plot(Xc, rho, Xv, Vx)
% Action
for it = 1:nt
    Vxc          = (Vx(2:end) + Vx(1:end-1))/2;
    rhoVx        = rho.*(Vx(2:end) + Vx(1:end-1))/2;
    rhoext       = [rho(1) (rho(2:end) + rho(1:end-1))/2 rho(end)] + 1e-8*Xv;
    Rrho         = (Vx0>0)*(rhoext(3:end)-rhoext(2:end-1))./(rhoext(2:end-1)-rhoext(1:end-2)) + (Vx0<0)*(rhoext(2:end-1)-rhoext(1:end-2))./(rhoext(3:end)-rhoext(2:end-1));
    Rrho(isnan(Rrho)==1) = 1;
    Rrho(Rrho==-Inf)     = -1;
    Rrho(Rrho==Inf)      = 1;
    phirho       = zeros(1,nx-1);
%     phirho(Rrho<1 & Rrho>0) = min(1,2*Rrho(Rrho<1 & Rrho>0)); 
%     phirho(Rrho>1)          = min(2,Rrho(Rrho>1));
%     phirho       = (Rrho+abs(Rrho))./(1+Rrho.^2);
    phirho(Rrho>0)          = min(1,Rrho(Rrho>0)); 
    %phirhoext    = 0.5*([phirho(1), phirho] + [phirho, phirho(end)]);
    Fcrho        = 0.5.*(rhoVx(2:end)+rhoVx(1:end-1)) - sign(Vx0)*0.5*diff(rhoVx) + 0.5*1.*(sign(Vx0)-Vx(2:end-1)*dt/dx).*diff(rhoVx) ;
    rho(2:end-1) = rho(2:end-1) - diff(Fcrho)/dx*dt;
    rho(end)     = rho(end) - (Vx0>0)*Vx0*(rho(end)-rho(end-1))/dx*dt;
    rho(1)       = rho(1)   - (Vx0<0)*Vx0*(rho(2)-rho(1))/dx*dt;
    rho2(2:end)  = rho2(2:end)   - (Vx0>0)*Vx0*diff(rho2)/dx*dt;
    rho2(1:end-1)= rho2(1:end-1) - (Vx0<0)*Vx0*diff(rho2)/dx*dt;
%     x2(2:end-1)   = Xc(2:end-1) - Vx0*dt;
%     rho2(2:end-1) = interp1( Xc,rho2,x2(2:end-1),'spline' );
    
    
    if mod(it-1, nout) == 0
%         plot(Xc, rho,'-k',Xc, rho2,'--r',Xc, rho0*heaviside(Xc-Vx0*dt*it)+1,'--b'), axis([-Lx/2 Lx/2 0 rho0+1.5 ]), drawnow
        plot(Xc, rho,'-x',Xc, rho2,'-r',Xc, rho0*exp(-3000*(Xc-Vx0*dt*it-0.3).^2)+1,'-b'), axis([-Lx/2 Lx/2 0 rho0+1.5 ]), drawnow
        %plot(phirho), drawnow
    end
end

[sum(rho*dx) sum(rho2*dx) sum((rho0*exp(-3000*(Xc-Vx0*dt*it-0.3).^2)+1)*dx)]-sum((rho0*exp(-3000*(Xc-Vx0*dt*it-0.3).^2)+1)*dx)

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


