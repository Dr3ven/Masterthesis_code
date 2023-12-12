clear,figure(1),clf % 1D shock wave
% Physics
Lx    = 1;
eta   = 0*1e-10;
V0    = 0;
rho0  = 1;
rhoa  = 10.05;
gam   = 1.4; 
lam   = Lx/10;
% Numerics
nx    = 5001;
nt    = 20*nx;
nout  = 1;
%preprocessing
dx    = Lx/(nx-1);
x     = -Lx/2:dx:Lx/2;
% Initial conditions
% rho   = rho0 + rhoa*exp(-(x/lam).^2);
rho(x>0) = rho0;  rho(x<=0) = rho0*3;
rhoVx = zeros(1,nx);
rhoE  = 1*rho;
rhoEVx= rhoVx;
Vx    = rhoVx./rho;
Pr    = (rhoE-Vx.^2/2)./(gam-1);
time  = 0;
w     = [rho;   rhoVx;          rhoE];
% Action
for it=1:600
    rho            = w(1,:);
    rhoVx          = w(2,:);
    rhoE           = w(3,:);
    Vx             = rhoVx./rho;
    Pr             = (rhoE-Vx.^2/2)./(gam-1);
    Vs             = sqrt(gam*Pr./rho) + 1*max(abs(Vx));
    dt             = min(dx^2/eta*max(rho)/4.1,dx/max(abs(Vs))/2.1);   
    F              = [rhoVx; Pr + rho.*Vx.^2; rhoE.*Vx + Vx.*Pr];
    Fc             = (F(:,2:end) + F(:,1:end-1))/2 ...
                   - (w(:,2:end) - w(:,1:end-1))/2*dx/dt;
    w(:,2:end-1)   = w(:,2:end-1) - diff(Fc,1,2)/dx*dt;
    time           = time + dt;
    if mod(it,nout)==0
        plot(x,rhoVx,'b',x,rho,'-rx',x,Pr,'k',x,rhoE,'g','Linewidth',1),title(it),drawnow
        drawnow
    end
end


