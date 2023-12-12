clear, figure(1), clf
% Physics
Lx      = 1e3;
rho0    = 2600;
beta    = 1e-11;

% numerics
nx      = 300;
dx      = Lx/(nx-1);
nout    = 10;
% Initialization
X             = -Lx/2:dx:Lx/2;
Xc            = -(Lx-dx)/2:dx:(Lx-dx)/2;
Xv            = -(Lx+dx)/2:dx:(Lx+dx)/2;

P       = 1e11*exp(-1e-3*X.^2);
%P       = 1e8*sin(X./(Lx/pi/10))./X;  
Pc      = 0.5*(P(2:end)+P(1:end-1));
rho     = rho0*(exp(beta.*Pc));
Vx      = zeros(1,nx+1);
P2      = P;
Vx2     = Vx;
rho2    = rho;

plot(P)
% Action
for it = 1:5e4
  dt = dx/sqrt(1/beta/min([rho(:); rho2(:)]))/6;
  
  dVxdt = -diff(P)/dx./rho;
  Vx(2:end-1) = Vx(2:end-1) + dVxdt*dt;
  Vx(1:end-1) = Vx(1:end-1) - (Vx(1:end-1)<0).*Vx(1:end-1).*diff(Vx)/dx*dt;
  Vx(2:end  ) = Vx(2:end  ) - (Vx(2:end  )>0).*Vx(2:end  ).*diff(Vx)/dx*dt;

  divV  = diff(Vx)/dx;
  dPdt  = -divV/beta;
  P     = P + dPdt*dt;
  Pc    = 0.5*(P(2:end)+P(1:end-1));
  rho   = rho0*(exp(beta.*Pc));
  
  dVxdt2 = -diff(P2)/dx./rho2; 
  Vx2(2:end-1) = Vx2(2:end-1) + dVxdt2*dt;
  divV2  = diff(Vx2)/dx;
  dPdt2  = -divV2/beta;
  P2     = P2 + dPdt2*dt;
  Pc2    = 0.5*(P2(2:end)+P2(1:end-1));
  rho2   = rho0*(exp(beta.*Pc2));

  
  if mod(it,nout) == 0
     figure(1), clf
     plot(Xc,Vx(2:end-1)./sqrt(1./(rho*beta)),Xc,Vx2(2:end-1)./sqrt(1./(rho2*beta))), axis 'tight', colorbar, drawnow
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


