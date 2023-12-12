clear,figure(1),clf % 1D shock wave
% Physics
Lx     = 1e4;       % Length, m 
E0     = 5e5;       % energy, J/kg  = Cp*T = 1000*500
rho0   = 3e3;       % density kg/m^3
rhoa   = 0.3*rho0;    % density kg/m^3
bet    = 3e-10;     % compressibility in 1/Pa
gru    = 1/30;      % alpha/beta/rho/Cp = 1/30, nondim
lam    = Lx/1000;   % initial shock width, m
Vximp  = 1e3;       % velocity of impact, m/s
% Numerics
nx     = 100;
nt     = 100000;
nout   = 1000;
%preprocessing
dx     = Lx/(nx-1);
x      = -Lx/2:dx:Lx/2;
dt     = dx/1e5;
% Initial conditions
%rho    = rho0 + rhoa*exp(-(x/lam).^2);
rho      = rho0 + rhoa*(1+tanh(-(x/lam)))/2;
w        = zeros(3,nx);
w(1,:)   = rho; 
w(2,x<0) = rho(x<0)*Vximp; Vx0 = w(2,:)./w(1,:);
w(3,:)   = rho.*(E0 + Vx0.^2/2); 
time     = 0;
% Action
for it=1:nt
    Vx             = w(2,:)./w(1,:);
    Pr             = (w(3,:) - w(2,:).*Vx/2)*gru  ...
                   - (1+gru)/bet*log(rho0./w(1,:));       
    F              = [w(2,:); Pr + w(2,:).^2./w(1,:); w(3,:).*Vx + Vx.*Pr];
    Fc             = (F(:,2:end) + F(:,1:end-1))/2 ...
                   - (w(:,2:end) - w(:,1:end-1))/2*dx/dt;
    w(:,2:end-1)   = w(:,2:end-1) - diff(Fc,1,2)/dx*dt;
    time           = time + dt;
    if mod(it,nout)==1
        subplot(311),plot(x,Vx0,x,Vx,'.'),ylabel('Vx')
        subplot(312),plot(x,Pr,'.'),ylabel('Pr')
        subplot(313),plot(x,w(1,:),'.'),ylabel('\rho')
        xlabel('x'),drawnow
    end
end


