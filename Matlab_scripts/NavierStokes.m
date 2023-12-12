clear, figure(1), clf
% Physics
Lx      = 1;
Ly      = 1;
R       = Lx/10;
eta     = 0.01;
rho     = 1;
beta    = 1;
dV      = 1e-2;
% numerics
nx      = 100;
ny      = 100;
dx      = Lx/(nx-1);
dy      = Ly/(ny-1);
nout    = 10;
% Initialization
X             = -Lx/2:dx:Lx/2;
Y             = -Ly/2:dy:Ly/2;
Xc            = -(Lx-dx)/2:dx:(Lx-dx)/2;
Yc            = -(Ly-dy)/2:dy:(Ly-dy)/2;
Xv            = -(Lx+dx)/2:dx:(Lx+dx)/2;
Yv            = -(Ly+dy)/2:dy:(Ly+dy)/2;
[x2d   y2d  ] = ndgrid(X ,Y );
[x2dc  y2dc ] = ndgrid(Xc,Yc);
[x2dVx y2dVx] = ndgrid(Xv,Y );
[x2dVy y2dVy] = ndgrid(X ,Yv);

P       = zeros(nx  ,ny  );
P       = exp(-100*(x2d.^2 + y2d.^2)); 
Vx      = zeros(nx+1,ny  );
Vy      = zeros(nx  ,ny+1);
Tauxy   = zeros(nx+1,ny+1);
dt = min(min(dx,dy)/sqrt(1/beta/max(rho(:)))/8.1);
imagesc(P)
% Action
for it = 1:5e4
  divV  = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
  dPdt  = -divV/beta;
  P     = P + dPdt*dt;
  dVxdt = -diff(P,1,1)/dx; 
  dVydt = -diff(P,1,2)/dy;
  Vx(2:end-1, :     ) = Vx(2:end-1, :     ) + dVxdt*dt;
  Vy( :     ,2:end-1) = Vy( :     ,2:end-1) + dVydt*dt;
  
  Vxc   = 0.5*(Vx(1:end-1,:)+Vx(2:end,:));
  Vyc   = 0.5*(Vy(:,1:end-1)+Vy(:,2:end));
  Vx(1:end-1,:) = Vx(1:end-1,:) - (Vx(1:end-1,:)<0).*Vx(1:end-1,:).*diff(Vx,1,1)/dx*dt;
  Vx(2:end  ,:) = Vx(2:end  ,:) - (Vx(2:end  ,:)>0).*Vx(2:end  ,:).*diff(Vx,1,1)/dx*dt;
%  Vx(2:end-1,1:end-1) = Vx(2:end-1,1:end-1) - (Vy(2:end-1,1:end-1)<0).*Vy(2:end-1,1:end-1).*diff(Vx(2:end-1,:),1,2)/dy*dt;
%  Vx(:,2:end  ) = Vx(:,2:end  ) - (Vy(:,2:end-1)>0).*Vy(:,2:end-1).*diff(Vxc,1,2)/dy*dt;
%   Vx(2:end  ,:) = Vx(2:end  ,:) - (Vx(2:end-1,:)>0).*Vx(2:end-1,:).*diff(Vxc,1,1)/dx*dt;        
%   Vx(:,1:end-1) = Vx(:,1:end-1) - (Vy(:,2:end-1)<0).*Vy(:,2:end-1).*diff(Vxc,1,2)/dy*dt;
%   Vx(:,2:end  ) = Vx(:,2:end  ) - (Vy(:,2:end-1)>0).*Vy(:,2:end-1).*diff(Vxc,1,2)/dy*dt;
%   
%   Vy(1:end-1,:) = Vx(1:end-1,:) - (Vx(2:end-1,:)<0).*Vx(2:end-1,:).*diff(Vxc,1,1)/dx*dt;
%   Vy(2:end  ,:) = Vx(2:end  ,:) - (Vx(2:end-1,:)>0).*Vx(2:end-1,:).*diff(Vxc,1,1)/dx*dt;        
%   Vy(:,1:end-1) = Vx(:,1:end-1) - (Vy(:,2:end-1)<0).*Vy(:,2:end-1).*diff(Vxc,1,2)/dy*dt;
%   Vy(:,2:end  ) = Vx(:,2:end  ) - (Vy(:,2:end-1)>0).*Vy(:,2:end-1).*diff(Vxc,1,2)/dy*dt;
  
  if mod(it,nout) == 0
     figure(1), clf
     surf(P), shading interp, axis 'tight', colorbar, drawnow
  end
  
end
figure(2)
subplot(221)
imagesc(flipud(P')), axis 'tight', colorbar, title('P'), drawnow
subplot(222)
imagesc(flipud(Vx')), axis 'tight', colorbar, title('Vx'), drawnow
subplot(223)
imagesc(flipud(Vy')), axis 'tight', colorbar, title('Vy'), drawnow

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


