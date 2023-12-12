clear, figure(2), clf
% Physics
Lx      = 1;
rho0    = 1;
Vx0     = 1e-2;
D       = 2e-5;
beta    = 1e-2;
mu      = 0*1e-3;
% numerics
nx      = 1501;
dx      = Lx/(nx);
nout    = 1;
c       = 1/sqrt(rho0*beta);
dtPT    = min([dx/abs(Vx0),dx/c])/3;
dt      = 1*dtPT;
CFL     = 8;
% Initialization
Xv      = -Lx/2:dx:Lx/2;
Xc      = -(Lx+dx)/2:dx:(Lx+dx)/2;

% rho     = 1.0005*exp(-300*(Xc+0).^2)+rho0;
% rho     = heaviside(Xc)+rho0;
rho(Xc>0) = rho0;  rho(Xc<=0) = rho0*15;
P       = 1/beta*log(rho(2:end-1)/rho0);
Vx      = 0*Vx0*ones(1,nx+1);
Mx      = av(rho).*Vx;
rhoc         = 0.5*(rho(1:end-1)+rho(2:end));
Vxc          = 0.5*(Vx(1:end-1)+Vx(2:end));
rho2 =rho; Vx2 =Vx; Mx2 = Mx; P2 =P; Vx2c = Vxc;
rho_ini =sum(rho(2:end-1)*dx);
Mx_ini = sum(Mx*dx);
P_ini = sum(P);
rhoRes = zeros(1,nx);
MxRes  = zeros(1,nx-1);
%Vx      = linspace(Vx0/2,Vx0,length(Xc));
plot(Xc, rho, Xc(2:end-1), P)
marker_pos = 0;
marker_neg = 0;
BForce = 0;
% Action
for it = 1:2e3

    rho_old = rho;
    Mx_old = Mx;
    P_old = P;
    rho2_old = rho2;
    Mx2_old = Mx2;
    Frhox2_old = (Vx2>0).*Vx2.*rho2(1:end-1) + (Vx2<0).*Vx2.*rho2(2:end);
    P2_old = 1/beta*log(rho2(2:end-1)/rho0);
    Sxx2_old = -P2 + mu.*diff(Vx2/dx);
    FMxx2_old = (Vx2c>0).*Vx2c.*Mx2(1:end-1) + (Vx2c<0).*Vx2c.*Mx2(2:end);
    err = 1;
    iter = 0;

    while err > 1e-8
        iter = iter+1;
        dt          = min([dx/max(abs(Vx)),dx/max(c)])*2;
        % BAckward Euler cons of mass
        Frhox        = (Vx>0).*Vx.*rho(1:end-1) + (Vx<0).*Vx.*rho(2:end);
%         Frhox        = Vx.*av(rho);
        drhodt       = (rho-rho_old)/dt;
        rhoRes_old   = rhoRes;
        rhoRes       = - drhodt(2:end-1) - diff(Frhox)/dx; 
        c_loc        = 1./sqrt(rho(2:end-1)*beta);
        dtPT_rho     = min(dx./abs(av(Vx))/CFL,dx./c_loc/CFL);
        rho(2:end-1) = rho(2:end-1) + rhoRes.*dtPT_rho;
        rho(1)       = rho(2);
        rho(end)     = rho(end-1);
        
        % Crank-Nicholson cons of mass     
        Frhox2        = (Vx2>0).*Vx2.*rho2(1:end-1) + (Vx2<0).*Vx2.*rho2(2:end);
        Frhox2        =  0.5*(Frhox2+Frhox2_old);
        drho2dt       = (rho2-rho2_old)/dt;
        rho2Res       = - drho2dt(2:end-1) - diff(Frhox2)/dx; 
        rho2(2:end-1) = rho2(2:end-1) + rho2Res*dtPT;
        rho2(1)       = rho2(2);
        rho2(end)     = rho2(end-1);

        % Backward Euler cons of momentum
        rhoc         = 0.5*(rho(1:end-1)+rho(2:end));
        Vxc          = 0.5*(Vx(1:end-1)+Vx(2:end));
        Mx           = rhoc.*Vx;
        % analytical pressure formula
        P            = 1/beta*log(rho(2:end-1)/rho0);
        %%%% incremental pressure formula %%%
%         dPdt         = (P-interp1(Xc(2:end-1),P_old,Xc(2:end-1)-dt*Vxc,'linear','extrap'))/dt;
%         dlnrhodt     = (log(rho(2:end-1))-interp1(Xc(2:end-1),log(rho_old(2:end-1)),Xc(2:end-1)-dt*Vxc,'linear','extrap'))/dt;
%         PRes         = - dPdt + 1/beta*dlnrhodt;
%         P            = P + PRes.*dtPT_rho;
        %%% incremental pressure formula %%%
        Sxx          = -P + mu.*diff(Vx)/dx;
%         FMxx         = (Vxc>0).*Vxc.*Mx(1:end-1) + (Vxc<0).*Vxc.*Mx(2:end);
        FMxx         = (Vxc>0).*Vxc.*Vxc.*rho(2:end-1) + (Vxc<0).*Vxc.*Vxc.*rho(2:end-1);
%         FMxx         =  Vxc.*av(Mx);
        dMxdt       = (Mx-Mx_old)/dt;
        MxRes_old   = MxRes;
        MxRes       = -dMxdt(2:end-1)-diff(FMxx-Sxx)/dx;
        dtPT_Mx     =   min(dx./abs(Vx(2:end-1))/CFL,av(dtPT_rho));
        Mx(2:end-1) = Mx(2:end-1) + MxRes.*dtPT_Mx;
%         Vxc_left    = -Vxc(1);
%         Mx_left     = -Mx(2);
%         FMxx_left   = (Vxc_left>0).*Vxc_left.*Mx_left;
%         P_left      = (FMxx(1)-Sxx(1))-FMxx_left;
%         Vxc_right    = -Vxc(end);
%         Mx_right     = -Mx(end-1);
%         FMxx_right   = (Vxc_right<0).*Vxc_right.*Mx_right;
%         P_right      = (FMxx(end)-Sxx(end))-FMxx_right;
% %         P_right     = -(FMxx(end)-Sxx(end));
%         rho(1)      = rho0*exp(beta*P_left);
%         rho(end)    = rho0*exp(beta*P_right);
% %         Mx(end)      = Mx(end) - (Vxc(end)>0)*Vxc(end)*(Mx(end)-Mx(end-1))/dx*dt;
% %         Mx(1)        = Mx(1  ) - (Vxc(1  )<0)*Vxc(1  )*(Mx(2  )-Mx(1    ))/dx*dt;
        Vx           = Mx./rhoc;
        % Crank-Nicholson cons of momentum
        rho2c         = 0.5*(rho2(1:end-1)+rho2(2:end));
        Vx2c          = 0.5*(Vx2(1:end-1)+Vx2(2:end));
        Mx2           = rho2c.*Vx2;
        P2            = 1/beta*log(rho2(2:end-1)/rho0);
        Sxx2          = -P2 + mu.*diff(Vx2/dx);
        FMxx2         = (Vx2c>0).*Vx2c.*Mx2(1:end-1) + (Vx2c<0).*Vx2c.*Mx2(2:end);
        FMxx2         =  0.5*(FMxx2+FMxx2_old);
        dMx2dt       = (Mx2-Mx2_old)/dt;
        Mx2Res       = -dMx2dt(2:end-1)-diff(FMxx2-Sxx2)/dx;

        Mx2(2:end-1)  = Mx2(2:end-1) + Mx2Res*dtPT;

    %     Mx(end)      = Mx(end) - (Vxc(end)>0)*Vxc(end)*(Mx(end)-Mx(end-1))/dx*dt;
    %     Mx(1)        = Mx(1  ) - (Vxc(1  )<0)*Vxc(1  )*(Mx(2  )-Mx(1    ))/dx*dt;
        Vx2           = Mx2./rho2c;
        
        if mod(iter-1, 25) == 0
            err = max(abs([rhoRes(:); MxRes(:)]));
%             semilogy(iter,err,'o'), hold on, drawnow
        end
    end
    marker_pos = marker_pos + 1/sqrt(rho0*beta)*dt;
    marker_neg = marker_neg - 1/sqrt(max(rho)*beta)*dt;
%     BForce = BForce + (1/beta*log((rho(1)+rho(2))/2/rho0)-1/beta*log((rho(end)+rho(end-1))/2/rho0))*dt;
    BForce = BForce  + ( (FMxx(1)-Sxx(1)) - (FMxx(end)-Sxx(end)) )*dt;
    if mod(it-1, nout) == 0
%         plot(Xc, rho,'-',Xv, Vx,'-',Xc, rho2,'-.',Xv, Vx2,'-.'), drawnow
        subplot(3,1,1:2)
        plot(Xc, rho,'-',Xc(2:end-1), P,'-',Xv, Vx,'-',Xv, Mx,'-',[marker_pos marker_pos],[-4 5],[marker_neg marker_neg],[-4 5])
        legend('rho','P','Vx','rhoVx')
        subplot(313)
        semilogy(it,abs(sum(rho(2:end-1)*dx)-rho_ini),'ob',it,abs(sum(Mx(2:end-1)*dx)-BForce),'or'),  drawnow, hold on
%         plot(Xv(2:end-1), (MxRes),'-'), drawnow
    end
end


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


