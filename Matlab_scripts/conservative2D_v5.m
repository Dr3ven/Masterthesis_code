clear, figure(1), %clf, colormap(parula(256))
% Physics
Lx      = 1e1;
Ly      = Lx/3;
rho0    = 1;
Vx0     = 1;
beta    = 1e-6;
mu      = 3e-3;
% characteristic values
% tc      = beta*mu*1e3;
% rhoc    = rho0;
% Vc      = Vx0;
% Lc      = Vc*tc;
% muc     = rhoc*Lc*Vc;
% betac   = 1/(muc/tc);
muc     = mu;
rhoc    = rho0;
betac   = beta;
tc      = muc*betac;
Lc      = sqrt(muc/rhoc*tc);
Vc      = Lc/tc;

% nondimensional values
Lx      = Lx/Lc;
Ly      = Ly/Lc;
rho0    = rho0/rhoc;
Vx0     = Vx0/Vc;
beta    = beta/betac;
mu      = mu/muc;
% numerics
nx      = 300;
ny      = 100;
dx      = Lx/(nx);
dy      = Ly/(ny);
nout    = 1;
dt      = dx/Vx0/2.1;
CFL_P   = 1;
CFL_V   = 1;
ksi     = 35;
% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx+dx)/2:dx:(Lx+dx)/2;
Yv      = -Ly/2:dy:Ly/2;
Yc      = -(Ly+dy)/2:dy:(Ly+dy)/2;
[x2dc  y2dc ] = ndgrid(Xc,Yc);
[x2dVx y2dVx] = ndgrid(Xv,Yc);
[x2dVy y2dVy] = ndgrid(Xc,Yv);

rho     = 0*exp(-300*(x2dc.^2 + y2dc.^2))+rho0;
% rho     = 0.25*heaviside(x2dc)*heaviside(y2dc)+rho0;
% rho     = 0.25*heaviside(x2dc)+rho0;
% rho(x2dc>0) = 1.25;  rho(x2dc<=0) = 0.75;
% rho(y2dc>0) = 1.25;  rho(y2dc<=0) = 0.75;
% rho = rho0*ones(size(x2dc)); rho(y2dc>0 & x2dc>0) = 16;
P         = 1/beta*log(rho(2:end-1)/rho0); 
Vx        = Vx0*ones(nx+1,ny+2);
Vx(1,:)   = Vx0;
Vx(end,:) = Vx0;
locX      = -Lx/2.5;
radrho    = sqrt((x2dc -locX).^2 + y2dc .^2);
radVx     = sqrt((x2dVx-locX).^2 + y2dVx.^2);
radVy     = sqrt((x2dVy-locX).^2 + y2dVy.^2);

diam      = Lx/10;
circ_X    = locX-diam/2:1e-4*Lx:locX+diam/2;
circ_Y    = sqrt((diam/2).^2 - (circ_X-locX).^2);
Vx(radVx<diam/2)       = 0;
Vy        = Vx0/1*exp(-1e3*((x2dVy-(locX+1.5*diam)   ).^2 + y2dVy.^2));
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

% hold off
% Action
for it = 1:1

    rho_old = rho;
    Mx_old = Mx;
    My_old = My;
    err = 1;
    iter = 0;
    dMxdtau   = 0*Mx(2:end-1,2:end-1);
    dMydtau   = 0*My(2:end-1,2:end-1);
    while err > 1e-8
        iter = iter+1;
%         dt           = min([dx/max(abs(Vx(:))),dy/max(abs(Vy(:))),min(dx,dy)*sqrt(max(rho(:))*beta),min(dx,dy).^2/mu])*2;
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
%         rhoRes(radrho(2:end-1,2:end-1)<0.1)       = 0;
        rho(2:end-1,2:end-1) = rho(2:end-1,2:end-1) + rhoRes.*dtrho.*CFL_P;
        % BC Inflow and outflow densities
        rho(1,:)       = rho(2,:);
        rho(end,:)     = rho(end-1,:);
        % BC impermeable walls
        rho(:,1)       = rho(:,2);
        rho(:,end)     = rho(:,end-1);
        
        % Strain-rates and stresses
        P            = 1/beta*log(rho/rho0);
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
        MxRes(radVx(2:end-1,2:end-1)<diam/2)       = 0;
        dMxdtau     = MxRes + 0*dMxdtau*(1-ksi/nx);
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
        MyRes       = - dMydt(2:end-1,2:end-1) - diff(FMyy - Syy,1,2)/dy - diff(FMyx - Sxy(:,2:end-1),1,1)/dx;
        MyRes(radVy(2:end-1,2:end-1)<diam/2)       = 0;
        dMydtau     = MyRes + 0*dMydtau*(1-ksi/ny);
        My(2:end-1,2:end-1)  = My(2:end-1,2:end-1) + dMydtau.*av_y(dtPT)*CFL_V;
        Vy           = My./av_y(rho);
        % BC fixed walls (normal velocity = 0)
        My(1,:)      = -My(2,:);
        My(end,:)    = -My(end-1,:);
        % BC no slip on horizontal walls
        My(:,1)      = My(:,2);
        My(:,end)    = My(:,end-1);

        if mod(iter, 1) == 0
            err = max(abs([rhoRes(:); MxRes(:); MyRes(:)]));
            if isnan(err) == 1 
                break
            end

%                 surf(x2dc,y2dc,rho), shading interp, drawnow
%                 imagesc(MxRes), colorbar, drawnow
            semilogy(iter,err,'o'), hold on, drawnow
        end
    end
    if mod(it-1, 1) == 0
%         plot(Xc, rho,'-',Xv, Vx,'-',Xc, rho2,'--',Xv, Vx2,'--'), axis([-Lx/2 Lx/2 -4*rho0 4*rho0]), drawnow
%         plot(Xc, rho,'-',Xv, Vx,'-'), axis([-Lx/2 Lx/2 -4*rho0 4*rho0]), drawnow
%         plot(it,sum(rho(2:end-1))-rho_ini,'ob',it,sum(Mx)-Mx_ini,'or',it,sum(rho2(2:end-1))-rho_ini,'xb',it,sum(Mx2)-Mx_ini,'xr'),  drawnow, hold on
%         plot(it,iter,'o'), hold on, drawnow
%     pcolor(x2dc,y2dc,rho), shading flat,title(it), hold on    
%     pcolor(x2dVy,y2dVy,Vy), shading flat,title(it), hold on
%     plot(circ_X,circ_Y,'k',circ_X,-circ_Y,'k')
%     X = av_xy(x2dc); Y = av_xy(y2dc); U= av_y(Vx); V=av_x(Vy);
%     step = 3;
%     quiver(X(1:step:end,1:step:end),Y(1:step:end,1:step:end),U(1:step:end,1:step:end),V(1:step:end,1:step:end),'w'), hold off, axis equal, drawnow
    
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


