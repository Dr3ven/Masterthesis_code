clear, figure(1), %clf, colormap(parula(256))
% Physics
Lx      = 1e4;
Ly      = Lx;
rho0    = 1.225;
drho    = 3e2;
Vx0     = 0;
P0      = 1e5;
beta    = 1/141e3;
mu      = 1.81e-5;
g       = 9.81;
% % characteristic values
% % tc      = beta*mu*1e3;
% % rhoc    = rho0;
% % Vc      = Vx0;
% % Lc      = Vc*tc;
% % muc     = rhoc*Lc*Vc;
% % betac   = 1/(muc/tc);
% 
% % muc     = mu;
% % rhoc    = rho0;
% % betac   = beta;
% % Pc      = 1/betac; 
% % tc      = muc*betac;
% % Lc      = sqrt(muc/rhoc*tc);
% % Vc      = Lc/tc;
% 
% Lc      = Lx;
% Pc      = P0;
% rhoc    = rho0+drho;
% betac   = 1/Pc;
% tc      = sqrt(rhoc*Lc^2/Pc);
% Vc      = Lc/tc;
% muc     = Pc*tc;
% 
% % Vc      = 1/sqrt(beta*rho0);
% % Lc      = Lx;
% % Pc      = P0;
% % tc      = Lc/Vc;
% % betac   = 1/Pc;
% % rhoc    = Pc/Vc^2;
% % 
% % muc     = Pc*tc;
% 
% % nondimensional values
% Lx      = Lx/Lc;
% Ly      = Ly/Lc;
% P0      = P0/Pc;
% rho0    = rho0/rhoc;
% drho    = drho/rhoc; 
% Vx0     = Vx0/Vc;
% beta    = beta/betac;
% mu      = mu/muc;
% numerics
nx      = 151;
ny      = 151;
dx      = Lx/(nx);
dy      = Ly/(ny);
nout    = 1;
dt      = 6e-2;
CFL_P   = 1/4;
CFL_V   = 1/16;
ksi     = 0*0.95;
% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx+dx)/2:dx:(Lx+dx)/2;
Yv      = -Ly/2:dy:Ly/2;
Yc      = -(Ly+dy)/2:dy:(Ly+dy)/2;
[x2dc  y2dc ] = ndgrid(Xc,Yc);
[x2dVx y2dVx] = ndgrid(Xv,Yc);
[x2dVy y2dVy] = ndgrid(Xc,Yv);

Vx        = Vx0*ones(nx+1,ny+2);
Vx(1,:)   = Vx0;
Vx(end,:) = Vx0;
locX      = 0;
locY      = -0.35*Ly;
diam      = 0.2*Lx;
radrho    = sqrt((x2dc -locX).^2 + (y2dc -locY).^2);
radVx     = sqrt((x2dVx-locX).^2 + (y2dVx-locY).^2);
radVy     = sqrt((x2dVy-locX).^2 + (y2dVy-locY).^2);
Xp        = [-1/2 -1/2 -0.3 -1/8 -0.02 -0.02  0.02  0.02  1/8  0.3  1/2  1/2]*Lx;
Yp        = [-1/2 -1/5 -1/5 -0.1 -0.15  -0.35 -0.35 -0.15  -0.1 -1/5 -1/5 -1/2]*Ly;
Xp2       = [-0.02 -0.02  0.02  0.02]*Lx;
Yp2       = [-0.1  -0.35 -0.35 -0.1]*Ly;
% plot(Xp,Yp)
maskrho  = 1-inpolygon(x2dc,y2dc,Xp,Yp);
maskrho(radrho<diam/2) = 1;
maskVx  = 1-inpolygon(x2dVx,y2dVx,Xp,Yp);
maskVx(radVx<diam/2) = 1;
maskVy  = 1-inpolygon(x2dVy,y2dVy,Xp,Yp);
maskVy(radVy<diam/2) = 1;

P         = P0*exp(-g*(y2dc+1/5*Ly)*0.028/288/8.314);
rho       = rho0.*exp(beta*(P-P0));
rho(radrho<diam/2) = rho0+drho;

rho(inpolygon(x2dc,y2dc,Xp2,Yp2)==1 & y2dc<locY+diam/2) = rho0+drho;
rho(inpolygon(x2dc,y2dc,Xp2,Yp2)==1 & y2dc>=locY+diam/2) = -(y2dc(inpolygon(x2dc,y2dc,Xp2,Yp2)==1 & y2dc>=locY+diam/2)+0.1*Ly)*drho/(-0.1-(locY+diam/2))+rho0; 
P         = 1/beta*log(rho/rho0)+P0; 
pcolor(x2dc,y2dc,P), shading flat,title('0'), colorbar, hold on
plot(polyshape(Xp,Yp))
Vy        = 0*Vx0/1*exp(-1e3*((x2dVy-(locX+1.5*diam)   ).^2 + y2dVy.^2));
Mx        = av_x(rho).*Vx;
My        = av_y(rho).*Vy;
rho_ini   = sum(rho(2:end-1));
Mx_ini    = sum(Mx);

rhoRes    = 0*rho(2:end-1,2:end-1);
%Vx      = linspace(Vx0/2,Vx0,length(Xc));
% pcolor(x2dVx,y2dVx,Vx), shading flat,title('0'), hold on
% plot(circ_X,circ_Y,'k',circ_X,-circ_Y,'k')
% X = av_xy(x2dc); Y = av_xy(y2dc); U= av_y(Vx); V=av_x(Vy);
% step = 3;
% quiver(X(1:step:end,1:step:end),Y(1:step:end,1:step:end),U(1:step:end,1:step:end),V(1:step:end,1:step:end),'w'), hold off, axis equal, drawnow

hold off
% Action
time = 0;
for it = 1:3e2
    rho_old = rho;
    Mx_old = Mx;
    My_old = My;
    err = 1;
    iter = 0;
    dMxdtau   = 0*Mx(2:end-1,2:end-1);
    dMydtau   = 0*My(2:end-1,2:end-1);
    while err > 1e-3
        iter = iter+1;
        dt           = min([dx/max(abs(Vx(:))),dy/max(abs(Vy(:))),min(dx,dy)*sqrt(max(rho(:))*beta),min(dx,dy).^2/mu])*8;
%         dt           = min([min(dx,dy)*sqrt(max(rho(:))*beta),min(dx,dy).^2/mu]);
        c_loc        = 1./sqrt(rho(2:end-1,2:end-1)*beta);
        dtPT         = min(min(min(dx./abs(av_x(Vx(:,2:end-1))),dx./abs(av_y(Vy(2:end-1,:)))),min(dx,dy).^2/mu*ones(nx,ny)),dx./c_loc);
        dtrho        = 1./(1/dt + 1./(min(dx,dy)./c_loc/4.1));
        % Conservation of mass
        Frhox        = (Vx>0).*Vx.*rho(1:end-1, :     ) + (Vx<0).*Vx.*rho(2:end, :   );
        Frhoy        = (Vy>0).*Vy.*rho( :     ,1:end-1) + (Vy<0).*Vy.*rho( :   ,2:end);
%         Frhox        = Vx.*av_x(rho);
%         Frhoy        = Vy.*av_y(rho);
        drhodt       = (rho-rho_old)/dt;
        rhoRes       = - drhodt(2:end-1,2:end-1) - diff(Frhox(:,2:end-1),1,1)/dx - diff(Frhoy(2:end-1,:),1,2)/dy; 
        rhoRes       = rhoRes.*maskrho(2:end-1,2:end-1);
        rho(2:end-1,2:end-1) = rho(2:end-1,2:end-1) + rhoRes.*dtrho.*CFL_P;
        % BC Inflow and outflow densities
        rho(1,:)       = rho(2,:);
        rho(end,:)     = rho(end-1,:);
        % BC impermeable walls
        rho(:,1)       = rho(:,2);
        rho(:,end)     = rho(:,end-1);
        
        % Strain-rates and stresses
        P            = 1/beta*log(rho/rho0)+P0;
        divV         = diff(Vx(:,2:end-1),1,1)/dx + diff(Vy(2:end-1,:),1,2)/dy;
        Exx          = diff(Vx(:,2:end-1),1,1)/dx - 1/3*divV;
        Eyy          = diff(Vy(2:end-1,:),1,2)/dx - 1/3*divV;
        Ezz          =                            - 1/3*divV;  
        Exy          = 0.5*(diff(Vy,1,1)/dx + diff(Vx,1,2)/dx);
        Sxx          = -P(2:end-1,2:end-1) + 2*mu.*Exx;
        Syy          = -P(2:end-1,2:end-1) + 2*mu.*Eyy;       
        Sxy          =                       2*mu.*Exy;
        Szz          = -P(2:end-1,2:end-1) + 2*mu.*Ezz; 
        dtV          = 1./(1/dt + 1./(min(dx,dy).^2./mu/4))*CFL_V;
        
        % Conservation of the x-component of momentum
        Mx           = av_x(rho).*Vx;
        FMxx         = (av_x(Vx( :     ,2:end-1))>0).*av_x(Vx( :     ,2:end-1)).*Mx(1:end-1,2:end-1) + (av_x(Vx( :     ,2:end-1))<0).*av_x(Vx( :     ,2:end-1)).*Mx(2:end  ,2:end-1);
%         FMxx         = av_x(Vx(:,2:end-1).*Mx(:,2:end-1));
        FMxy         = (av_x(Vy(2:end-1, :     ))>0).*av_x(Vy(2:end-1, :     )).*Mx(2:end-1,1:end-1) + (av_x(Vy(2:end-1, :     ))<0).*av_x(Vy(2:end-1, :     )).*Mx(2:end-1,2:end  );

        dMxdt       = (Mx-Mx_old)/dt;
        MxRes       = - dMxdt(2:end-1,2:end-1) - diff((FMxx - Sxx),1,1)/dx - diff(FMxy - Sxy(2:end-1,:),1,2)/dy;
        MxRes       = MxRes.*maskVx(2:end-1,2:end-1);
        dMxdtau     = MxRes + dMxdtau*ksi;
        Mx(2:end-1,2:end-1)  = Mx(2:end-1,2:end-1) + dMxdtau.*av_x(dtPT)*CFL_V;
        % BC fixed walls (normal velocity = 0)
        Mx(:,1)      = Mx(:,2);
        Mx(:,end)    = Mx(:,end-1);
        % BC no slip on vertical walls
%         Mx(1,:)      = Mx(2,:);
%         Mx(end,:)    = Mx(end-1,:);
        
        Vx           = Mx./av_x(rho);
        Vx(1,:)      = Vx0;
        Vx(end,:)    = Vx0;
        
        % Conservation of the y component of momentum
        My           = av_y(rho).*Vy;
        FMyy         = (av_y(Vy(2:end-1, :     ))>0).*av_y(Vy(2:end-1, :     )).*My(2:end-1,1:end-1) + (av_y(Vy(2:end-1, :     ))<0).*av_y(Vy(2:end-1, :     )).*My(2:end-1,2:end  );
        FMyx         = (av_y(Vx( :     ,2:end-1))>0).*av_y(Vx( :     ,2:end-1)).*My(1:end-1,2:end-1) + (av_y(Vx( :     ,2:end-1))<0).*av_y(Vx( :     ,2:end-1)).*My(2:end  ,2:end-1);

        dMydt       = (My-My_old)/dt;
        MyRes       = - dMydt(2:end-1,2:end-1) - diff(FMyy - Syy,1,2)/dy - diff(FMyx - Sxy(:,2:end-1),1,1)/dx - g*av_y(rho(2:end-1,2:end-1));
        MyRes       = MyRes.*maskVy(2:end-1,2:end-1);
        dMydtau     = MyRes + dMydtau*ksi;
        My(2:end-1,2:end-1)  = My(2:end-1,2:end-1) + dMydtau.*av_y(dtPT)*CFL_V;
        Vy           = My./av_y(rho);
        % BC fixed walls (normal velocity = 0)
        My(1,:)      = -My(2,:);
        My(end,:)    = -My(end-1,:);
        % BC no slip on horizontal walls
        My(:,1)      = My(:,2);
        My(:,end)    = My(:,end-1);

        if mod(iter, 25) == 0
            err = max(abs([rhoRes(:); MxRes(:); MyRes(:)]));
            if isnan(err) == 1 
                break
            end

%                 surf(x2dc,y2dc,rho), shading interp, drawnow
%                 imagesc(MxRes), colorbar, drawnow
%             semilogy(iter,err,'o'), hold on, drawnow
        end
    end
    time = time +dt;
    if mod(it-1, 1) == 0
%         plot(Xc, rho,'-',Xv, Vx,'-',Xc, rho2,'--',Xv, Vx2,'--'), axis([-Lx/2 Lx/2 -4*rho0 4*rho0]), drawnow
%         plot(Xc, rho,'-',Xv, Vx,'-'), axis([-Lx/2 Lx/2 -4*rho0 4*rho0]), drawnow
%         plot(it,sum(rho(2:end-1))-rho_ini,'ob',it,sum(Mx)-Mx_ini,'or',it,sum(rho2(2:end-1))-rho_ini,'xb',it,sum(Mx2)-Mx_ini,'xr'),  drawnow, hold on
%         plot(it,iter,'o'), hold on, drawnow
    pcolor(x2dc,y2dc,P), shading flat,title(time), caxis([P0 P0*2]), colorbar, hold on 
    plot(Xp,Yp,'w')
%     pcolor(x2dVy,y2dVy,Vy), shading flat,title(it), hold on
%     plot(circ_X,circ_Y,'k',circ_X,-circ_Y,'k')
    X = av_xy(x2dc); Y = av_xy(y2dc); U= av_y(Vx); V=av_x(Vy);
    step = 5;
    quiver(X(1:step:end,1:step:end),Y(1:step:end,1:step:end),U(1:step:end,1:step:end),V(1:step:end,1:step:end),'w'), hold off, axis equal, drawnow
    
    end
end


function A = av_x(B)
    A = 0.5*(B(2:end,:)+B(1:end-1,:));
end
function A = av_y(B)
    A = 0.5*(B(:,2:end)+B(:,1:end-1));
end
function A = av_xy(B)
    A = 0.25*(B(2:end,2:end)+B(1:end-1,2:end)+B(2:end,1:end-1)+B(1:end-1,1:end-1));
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


