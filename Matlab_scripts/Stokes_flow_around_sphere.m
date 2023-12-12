clear, figure(1), clf
% Physics
Lx      = 1;
Ly      = 1;
R       = Lx/10;
mu_mat  = 0.1;
mu_inc  = mu_mat/1000;
rho_mat = 1;
rho_inc = 0.9;
gx      = 0;
gy      = 0;
beta    = 1/(30*mu_mat^2);
Vbc     = 1;
% numerics
nx      = 50;
ny      = 50;
dx      = Lx/(nx-1);
dy      = Ly/(ny-1);
nout    = 1000;
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
mu      = mu_mat*ones(size(x2d));
muc     = mu_mat*ones(size(x2dc));
rho     = rho_mat*ones(size(x2d));
rho(max(abs(x2d),abs(y2d))<(0.1*Lx))     =  rho_inc;
P       = rho.*gy.*y2d;
Vx      = zeros(nx+1,ny  );
Vy      = ones(nx  ,ny+1);
Tauxy   = zeros(nx+1,ny+1);
dt = min(min(dx,dy)/sqrt(1/beta/max(rho(:)))/8.1, min(dx,dy)^2/max(mu(:)/max(rho(:)))/8.1);

% Action
  divV  = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
for iter = 1:5e5
  dPdt  = -divV/beta;
  P     = P + dPdt*dt; 
  Tauxx = 2*mu.*(diff(Vx,1,1)/dx-divV/3);
  Tauyy = 2*mu.*(diff(Vy,1,2)/dy-divV/3);
  Tauxy(2:end-1,2:end-1) = muc.*(diff(Vx(2:end-1,:),1,2)/dy + diff(Vy(:,2:end-1),1,1)/dx);
  dVxdt = (-diff(P,1,1)/dx + diff(Tauxx,1,1)/dx + diff(Tauxy(2:end-1,:),1,2)/dy) + 0.5*(rho(1:end-1,:) + rho(2:end,:)).*gx; 
  dVydt = (-diff(P,1,2)/dy + diff(Tauyy,1,2)/dy + diff(Tauxy(:,2:end-1),1,1)/dx) + 0.5*(rho(:,1:end-1) + rho(:,2:end)).*gy;
  Vx(2:end-1, :     ) = Vx(2:end-1, :     ) + dVxdt*dt;
  Vy( :     ,2:end-1) = Vy( :     ,2:end-1) + dVydt*dt;
  Vx(max(abs(x2dVx),abs(y2dVx))<(0.1*Lx))     =  0;
  Vy(max(abs(x2dVy),abs(y2dVy))<(0.1*Ly))     =  0;
    divV  = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
  if mod(iter,nout) == 0
     figure(1)
     semilogy(iter,max(abs(divV(:))),'rx',iter,max(abs(dVxdt(:))),'k+',iter,max(abs(dVydt(:))),'ob'), hold on, drawnow
    %figure(2), clf
    %imagesc(flipud(dVxdt')), axis 'tight', colorbar, drawnow

  end
end
figure(2)
subplot(221)
imagesc(flipud(P')), axis 'tight', colorbar, title('P'), drawnow
subplot(222)
imagesc(flipud(Vx')), axis 'tight', colorbar, title('Vx'), drawnow
subplot(223)
imagesc(flipud(Vy')), axis 'tight', colorbar, title('Vy'), drawnow
subplot(224)
imagesc(flipud(rho')), axis 'tight', colorbar, title('rho'), drawnow, hold on
%subplot(224)
%imagesc(x2d,y2d,flipud(P')), axis 'tight', colorbar, title(iter), drawnow

























