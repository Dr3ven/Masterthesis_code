clear, figure(1), %clf, colormap(parula(256))
% Physics
Lx      = 1;
Ly      = 1;
rho0    = 1;
Vx0     = 1e-2;
beta    = 1e-1;
mu      = 1e-2;
% numerics
nx      = 151;
ny      = nx;
dx      = Lx/(nx);
dy      = Ly/(ny);
nout    = 1;
c       = 1/sqrt(rho0*beta);
dt      = min([dx/abs(Vx0),dx/c])/2;
CFL     = 16;
% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx+dx)/2:dx:(Lx+dx)/2;
Yv      = -Ly/2:dy:Ly/2;
Yc      = -(Ly+dy)/2:dy:(Ly+dy)/2;
[x2dc  y2dc ] = ndgrid(Xc,Yc);
[x2dVx y2dVx] = ndgrid(Xc,Yv);
[x2dVy y2dVy] = ndgrid(Xv,Yc);

% rho     = 8*exp(-300*(x2dc.^2 + y2dc.^2))+rho0;
% rho     = 0.25*heaviside(x2dc)*heaviside(y2dc)+rho0;
% rho     = 0.25*heaviside(x2dc)+rho0;
% rho(x2dc>0) = 1.25;  rho(x2dc<=0) = 0.75;
% rho(y2dc>0) = 1.25;  rho(y2dc<=0) = 0.75;
rho = rho0*ones(size(x2dc)); rho(y2dc>0 & x2dc>0) = 16;
P       = 1/beta*log(rho(2:end-1)/rho0);
Vx      = 0*Vx0*ones(nx+1,ny+2);
Vy      = 0*Vx0*ones(nx+2,ny+1);
Mx      = av_x(rho).*Vx;
My      = av_y(rho).*Vy;
rho_ini =sum(rho(2:end-1));
Mx_ini = sum(Mx);
%Vx      = linspace(Vx0/2,Vx0,length(Xc));
surf(x2dc,y2dc,rho), shading flat, colorbar
hold off
% Action
for it = 1:100

    rho_old = rho;
    Mx_old = Mx;
    My_old = My;
    err = 1;
    iter = 0;

    while err > 1e-6
        iter = iter+1;
        dt           = min([dx/max(abs(Vx(:))),dx/max(abs(Vy(:))),dx/sqrt(max(rho(:))*beta)]);
        c_loc        = 1./sqrt(rho(2:end-1,2:end-1)*beta);
        dtPT         = min(min(dx./abs(av_x(Vx(:,2:end-1))),dx./abs(av_y(Vy(2:end-1,:)))),dx./c_loc)/CFL;
%         dtPT         = dt*dtPT./dtPT/8;
        % Conservation of mass
        Frhox        = (Vx>0).*Vx.*rho(1:end-1, :     ) + (Vx<0).*Vx.*rho(2:end, :   );
        Frhoy        = (Vy>0).*Vy.*rho( :     ,1:end-1) + (Vy<0).*Vy.*rho( :   ,2:end);
%         Frhox        = Vx.*av_x(rho);
%         Frhoy        = Vy.*av_y(rho);
        drhodt       = (rho-rho_old)/dt;
        rhoRes       = - drhodt(2:end-1,2:end-1) - diff(Frhox(:,2:end-1),1,1)/dx - diff(Frhoy(2:end-1,:),1,2)/dy; 

        rho(2:end-1,2:end-1) = rho(2:end-1,2:end-1) + rhoRes.*dtPT;
        % BC impermeable walls
        rho(1,:)       = rho(2,:);
        rho(end,:)     = rho(end-1,:);
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

        % Conservation of the x-component of momentum
        Mx           = av_x(rho).*Vx;
        FMxx         = (av_x(Vx( :     ,2:end-1))>0).*av_x(Vx( :     ,2:end-1)).*Mx(1:end-1,2:end-1) + (av_x(Vx( :     ,2:end-1))<0).*av_x(Vx( :     ,2:end-1)).*Mx(2:end  ,2:end-1);
%         FMxx         = av_x(Vx(:,2:end-1).*Mx(:,2:end-1));
        FMxy         = (av_x(Vy(2:end-1, :     ))>0).*av_x(Vy(2:end-1, :     )).*Mx(2:end-1,1:end-1) + (av_x(Vy(2:end-1, :     ))<0).*av_x(Vy(2:end-1, :     )).*Mx(2:end-1,2:end  );

        dMxdt       = (Mx-Mx_old)/dt;
        MxRes       = - dMxdt(2:end-1,2:end-1) - diff((FMxx - Sxx),1,1)/dx - diff(FMxy - Sxy(2:end-1,:),1,2)/dy;
        Mx(2:end-1,2:end-1)  = Mx(2:end-1,2:end-1) + MxRes.*av_x(dtPT);
        % BC fixed walls (normal velocity = 0)
%         Mx(:,1)      = -Mx(:,2);
%         Mx(:,end)    = -Mx(:,end-1);
        % BC no slip on vertical walls
%         Mx(1,:)      = -Mx(2,:);
%         Mx(end,:)    = -Mx(end-1,:);
        Vx           = Mx./av_x(rho);

        
        % Conservation of the y component of momentum
        My           = av_y(rho).*Vy;
        FMyy         = (av_y(Vy(2:end-1, :     ))>0).*av_y(Vy(2:end-1, :     )).*My(2:end-1,1:end-1) + (av_y(Vy(2:end-1, :     ))<0).*av_y(Vy(2:end-1, :     )).*My(2:end-1,2:end  );
        FMyx         = (av_y(Vx( :     ,2:end-1))>0).*av_y(Vx( :     ,2:end-1)).*My(1:end-1,2:end-1) + (av_y(Vx( :     ,2:end-1))<0).*av_y(Vx( :     ,2:end-1)).*My(2:end  ,2:end-1);

        dMydt       = (My-My_old)/dt;
        MyRes       = - dMydt(2:end-1,2:end-1) - diff(FMyy - Syy,1,2)/dy - diff(FMyx - Sxy(:,2:end-1),1,1)/dx;
%         MyRes       = MyRes*0;
        My(2:end-1,2:end-1)  = My(2:end-1,2:end-1) + MyRes.*av_y(dtPT);
        Vy           = My./av_y(rho);
        % BC fixed walls (normal velocity = 0)
%         My(1,:)      = -My(2,:);
%         My(end,:)    = -My(end-1,:);
        % BC no slip on horizontal walls
%         My(:,1)      = -My(:,2);
%         My(:,end)    = -My(:,end-1);
        if mod(iter-1, 100) == 0
            err = max(abs([rhoRes(:); MxRes(:); MyRes(:)]));
%                 surf(x2dc,y2dc,rho), shading interp, drawnow
%                 imagesc(rho), drawnow
%             semilogy(iter,err,'o'), hold on, drawnow
        end
    end
    if mod(it-1, nout) == 0
%         plot(Xc, rho,'-',Xv, Vx,'-',Xc, rho2,'--',Xv, Vx2,'--'), axis([-Lx/2 Lx/2 -4*rho0 4*rho0]), drawnow
%         plot(Xc, rho,'-',Xv, Vx,'-'), axis([-Lx/2 Lx/2 -4*rho0 4*rho0]), drawnow
%         plot(it,sum(rho(2:end-1))-rho_ini,'ob',it,sum(Mx)-Mx_ini,'or',it,sum(rho2(2:end-1))-rho_ini,'xb',it,sum(Mx2)-Mx_ini,'xr'),  drawnow, hold on
%         plot(it,iter,'o'), hold on, drawnow
    surf(x2dc,y2dc,rho), shading interp, drawnow
%     imagesc(rho), drawnow
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


