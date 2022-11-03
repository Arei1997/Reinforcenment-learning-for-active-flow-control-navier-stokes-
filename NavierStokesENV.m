clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set:
% (0) decide whether restar (restart=1) or scratch (restart=0)
% (1) the number of points (number of pressure cells)
% (2) the dimension of the domain
% (3) the Reynolds number and other geometrical characteristics
% (4) (for the moment) a fixed time step
%      and the number of time steps
% (5) the boundary values and rkpm yes/no
% (6) check convergence every iconv time steps
%     eventually save the field and
%     decide whether or not open "on the fly" figures      
%                                                                   
%               Ub(3) and Vb(3)                                      
%           ¡---------------------¡                                   
%   Ub(4)   ¡                     ¡ Ub(2) and Vb(2)                    
%   Vb(4)   ¡                     ¡                                     
%           -----------------------                                      
%               Ub(1) and Vb(1)                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (0)
restart=0 % decide if from scratch (0) or from former solution (1)
% (1)
ncx=500; % number of nodes for the pressure in x
ncy=420; % number of nodes for the pressure in y
% (2)
lx=16.5; % x-wise extension of the domain 
ly=13.86; % y-wise extension of the domain 
% (3)
%Re=100; % Reynolds number of the immersed circle (U upstream * Diameter/kinematic visc.)
%Re=90;
%Re=80;
%Re=70;
%Re=60;
%Re=50;
%Re=130;
Re=150;
D=1;   %Diameter of the cylinder;
rc=D/2;% circle radius
xc=3.003;  % center of circle (coordinates)
yc=6.930;
% (4)
nt=1000; CFL=0.2; % number of time steps and CFL number
[Ub,Vb]=set_velocity_bc(ncx,ncy,lx,ly); % set velocity on the boundaries
% decide if rkpm or not
irkpm=0 % no rkpm
%irkpm=1 % yes rkpm
%irkpm_press=1 % (6) rkpm for pressure too 
irkpm_press=0 % (6) rkpm for pressure too 
iflap=1;
%iflap=1; % flapping rod
restart_flap=0;
imovie=0;
%imovie=1;
% set the initial flap position and constants
length_beam=D; center_rot=length_beam; 
k_spring=0.5; rho_beam=2; Delta_rho=rho_beam-1;
Inertia_mom=Delta_rho*length_beam^3/3;
Inertia_mom=Inertia_mom+Delta_rho*length_beam*center_rot*(center_rot-length_beam);
%%%% rotation coeficients
coef_rot1=k_spring/Inertia_mom;
coef_rot2=1/Inertia_mom;
%
iconv=20; % every iconv steps gives some info on the screen
snapshot_it=100; % every snapshot save the current solution
nsteps=[1:snapshot_it:nt];
yes_snap=0; % (0) no plots on the fly (1) yes
% assemble the grid and the initial solution
angolo=pi/100; angolo_dot=0; TORQUE=0;
%angolo=-pi/6; angolo_dot=0; TORQUE=0;
if restart==0
   [xu,xv,xp,yu,yv,yp,u,v,p]=initialize(ncx,ncy,lx,ly,Ub,Vb); %initialize the field and grid if from scratch
   save xp.mat xp;
   save yp.mat yp;
   save xu.mat xu;
   save yu.mat yu;
   save xv.mat xv;
   save yv.mat yv;
else % otherwise load previous solution and grid
   load uu; load vv; load pp;
   load xp; load yp; load xu; load xv; load yu; load yv;
   if restart_flap==1
      load angolo.mat angolo; load angolo_dot.mat angolo_dot; load TORQUE.mat TORQUE;
   end
end
% initialize the matrices for Helmholz and Poisson Problems
nxu=length(xu); nyu=length(yu); %no. of nodes for u (x-comp. of velocity) % CAREFUL STAGGERED GRID
nxv=length(xv); nyv=length(yv); %no. of nodes for v (y-comp. of velocity)
nxp=length(xp); nyp=length(yp); %no. of nodes for p 
% no uniform mesh with lx assumed not to be equal to ly!!!!
hx=lx/ncx; hy=ly/ncy; %mesh spacing
Deltat=CFL*hx/max(Ub(4,:)); % Allowed time step using entrance velocity as max
disp('Delta t ')
disp(Deltat)
beta=Deltat/hx; % some parameters
alpha=Re/Deltat; 
%h=xp(2)-xp(1);
% prefactor some parts of approximate factorization of viscous terms 
% We use D'Yakonov scheme (approximate factorization)
tic
[alxu,bexu,cxu,alyu,beyu,cyu]=predyakonovU(nxu,nyu,hx,hy,Deltat,Re);
[alxv,bexv,cxv,alyv,beyv,cyv]=predyakonovV(nxv,nyv,hx,hy,Deltat,Re);
[PX,PY,PMIX]=assmatFD(ncx,ncy,hx,hy); % diagonalize pressure laplacian operator
toc
disp('finished eig problems ')
%[xlag,ylag,theta,dS]=initialize_lag(xu,yu,rc,xc,yc,irkpm,0); % discretize circle with uniformly spaced Lagrange nodes
[xlag,ylag,theta,dS,nlsts,nlste]=initialize_lag_tail(xu,yu,rc,xc,yc,irkpm,0,angolo); 
%[xlag,ylag,dS,nlsts,nlste]=initialize_lag_onlytail(xu,yu,rc,xc,yc,center_rot,angolo); 
% compute the delta mollifier and the width of the corona around the cirlce (\epsilon in Jou. Comp. Phys. Paper)
if irkpm > 0
tic 
[epsilonx,bcoefx,ix_x,jy_x,area_x,nsup_x,hx_x,hy_x]=precompute_eps(xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,dS);
%keyboard
[epsilony,bcoefy,ix_y,jy_y,area_y,nsup_y,hx_y,hy_y]=precompute_eps(xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,dS);
toc
if irkpm_press==1
[epsilonp,bcoefp,ix_p,jy_p,area_p,nsup_p,hx_p,hy_p]=precompute_eps(xp(2:length(xp)-1),yp(2:length(yp)-1),xlag,ylag,dS);
end
else
[epsilonx,bcoefx,ix_x,jy_x,area_x,nsup_x,hx_x,hy_x]=precompute_eps_norkpm(xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,dS);
[epsilony,bcoefy,ix_y,jy_y,area_y,nsup_y,hx_y,hy_y]=precompute_eps_norkpm(xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,dS);
end
disp('finished precomputing immersed boundary data')
%keyboard
%%%% start the time loop
disp('starting time advancement')
isnapshot=0;
Qin=sum(u(nxu,2:nyu-1)); % flow rate at inlet
AB1=1; AB2=0; CGP=0; % first step Euler for nonlinear terms, from the second one Adams-Bashfort
laplax=zeros(nxu-2,nyu-2); NLx=zeros(nxu-2,nyu-2);
laplay=zeros(nxv-2,nyv-2); NLy=zeros(nxv-2,nyv-2);
icount=0;
plot(xlag,ylag,'k',xc,yc,'r*'); 
axis([xc-7*D/6 xc+2.5*D yc-D/2 yc+D/2]); hold on;
axis equal
wkeyboard; 
for it=1:nt
disp('step number'); disp(it)
tic
    actime=(it-1)*Deltat;
    if mod(it,iconv) == 0
       uold=u; vold=v;
    end
    if mod(it,10*iconv) == 0
      AB1=1; AB2=0;
    end

%%%%%%%% X_MOM EQUATION %%%%%%%%%%%%%%
% FROM A to B prediction of u from x momentum equation without immersed object (fully explicit)
% --------- A --------------
% compute the laplacian of u
   lapla=laplacian(u,hx,hy);
% compute NL terms of X-mom equation and dp/dx
   NL=nonlinear(u,v,Ub,Vb,lx,ly,ncx,ncy,1);
   GP=pressgrad(p,hx,hy,1);
% compute the force for the x-mom equation
   rhs=sumupAB(u,lapla,NL,laplax,NLx,CGP*GP,Deltat,Re,AB1,AB2); %(2.10a)
   if iflap==0
   [rhs,val_x]=immersed_boundary(rhs,xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,dS,epsilonx,bcoefx,...
       ix_x,jy_x,area_x,nsup_x,hx_x,hy_x); 
   else
   %velx=-angolo_dot*(ylag(nlsts:nlste)-ylag(nlsts));
   velx=-angolo_dot*(ylag(nlsts:nlste)-yc);
   [rhs,val_x]=immersed_boundary(rhs,xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,dS,epsilonx,bcoefx,...
       ix_x,jy_x,area_x,nsup_x,hx_x,hy_x,nlsts,nlste,velx/Deltat); 
   end
   FX=-val_x(nlsts:nlste).*epsilonx(nlsts:nlste).*dS(nlsts:nlste)';
% --------- B -------------- 
% rhs contains the spread force field required to impose bc's on the circle for u
%keyboard
   forcex=hx*hy*sum(sum(rhs)); % compute integral force in x direction on the circle
   Cd=2*forcex/D; %...and the drag coeficient
% assemble rhs for the implicit step
% from C to D redo x-momentum with the spreaded body force on the rhs 
% Approximate fact. for viscous terms and Adams Bashfort for the non-linear ones
%---------C ------------
   rhs=rhs-AB1*NL+AB2*NLx-CGP*GP;
   [Ubs]=setBCout(u,v,Deltat,hx,Qin,1);
   ustar=DyakonovU(rhs,u,Ub,Ubs,alxu,bexu,cxu,alyu,beyu,cyu,hx,hy,Re,Deltat);
%---------D ------------
   laplax=lapla; NLx=NL;
% DO STAGE A-B and C-D for the y momentum equation
%%%%%%%% Y_MOM EQUATION %%%%%%%%%%%%%%
% compute the laplacian of v
   lapla=laplacian(v,hx,hy);
% compute NL terms of Y-mom equation and dp/dy
   NL=nonlinear(u,v,Ub,Vb,lx,ly,ncx,ncy,2);
   GP=pressgrad(p,hx,hy,2);
% compute the force for the y-mom equation
   rhs=sumupAB(v,lapla,NL,laplay,NLy,CGP*GP,Deltat,Re,AB1,AB2);
   if iflap==0
   [rhs,val_y]=immersed_boundary(rhs,xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,dS,epsilony,bcoefy,...
             ix_y,jy_y,area_y,nsup_y,hx_y,hy_y); 
   else
   %vely=angolo_dot*(xlag(nlsts:nlste)-xlag(nlsts));
   vely=angolo_dot*(xlag(nlsts:nlste)-xc);
   [rhs,val_y]=immersed_boundary(rhs,xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,dS,epsilony,bcoefy,...
             ix_y,jy_y,area_y,nsup_y,hx_y,hy_y,nlsts,nlste,vely/Deltat); 
   end
   FY=-val_y(nlsts:nlste).*epsilony(nlsts:nlste).*dS(nlsts:nlste)';
    if mod(it,iconv) == 0
       disp('saving force field')
       save rhs.mat rhs;
    end
   forcey=hx*hy*sum(sum(rhs));
   Cf=2*forcey/D;
% assemble rhs for the implicit step
   rhs=rhs-AB1*NL+AB2*NLy-CGP*GP;
   [Vbs]=setBCout(u,v,Deltat,hx,Qin,2);
%   vstar=DyakonovV(rhs,v,Vb,alxv,bexv,cxv,alyv,beyv,cyv,hx,hy,Re,Deltat);
   vstar=DyakonovV(rhs,v,Vb,Vbs,alxv,bexv,cxv,alyv,beyv,cyv,hx,hy,Re,Deltat);
   laplay=lapla; NLy=NL;
%%%%%% Get ready for the projection step
% set the boundary values on the boundary frame
   [u,v]=setbc(ustar,vstar,Ub,Vb,Ubs,Vbs);
% prepare the rhs for the pressure equation
    div=divergence(u,v,hx,hy,ncx,ncy);
% solve the Poisson problem for pressure
   phi=solveFD(PX,PY,PMIX,div/Deltat);
   [u,v,p]=finalvelocityVK(phi,u,v,p,Deltat,Re,hx,hy);
   ulag=interpolante(u(2:nxu-1,2:nyu-1),xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,...
                         dS,bcoefx,ix_x,jy_x,area_x,nsup_x,hx_x,hy_x);
   vlag=interpolante(v(2:nxv-1,2:nyv-1),xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,...
                          dS,bcoefy,ix_y,jy_y,area_y,nsup_y,hx_y,hy_y);
v1(1)=max(abs(ulag));
v1(2)=max(abs(vlag));
fprintf(1,'Max err U=%8.3e ; Max err V (NC)=%8.3e \n',v1(1),v1(2));
  if mod(it,iconv) == 0
     icount=icount+1;
     disp('saving the field')
     isnapshot=isnapshot+1;
     dummy=screenout(uold,u,vold,v,actime,hx,hy,ncx,ncy,Cd,Cf);
     forcexx(isnapshot)=forcex; forceyy(isnapshot)=forcey;
     Cdd(isnapshot)=Cd; Cll(isnapshot)=Cf;
     time(isnapshot)=actime;
     dummy=savefield(u,v,p,phi,time,forcexx,forceyy,Cdd,Cll);
     if iflap==1
        ango(isnapshot)=angolo;
        ango_dot(isnapshot)=angolo_dot;
        torques(isnapshot)=TORQUE;
        save angolo.mat angolo; save angolo_dot.mat angolo_dot; save TORQUE.mat TORQUE;
        save ango.mat ango; save ango_dot.mat ango_dot; save torques.mat torques;
     end
     disp('Force X'); disp(sum(FX));
     disp('Force Y'); disp(sum(FY));
%     disp('angle'); disp(angolo);
%     disp('angular velocity'); disp(angolo_dot);
     if irkpm_press==1
     dummy=savelagran(u(2:nxu-1,2:nyu-1),v(2:nxv-1,2:nyv-1),div(2:nxp-1,2:nyp-1),...
           xu(2:nxu-1),yu(2:nyu-1),xv(2:nxv-1),yv(2:nyv-1),xp(2:nxp-1),yp(2:nyp-1),...
           xlag,ylag,dS,bcoefx,bcoefy,bcoefp,ix_x,ix_y,ix_p,jy_x,jy_y,jy_p,area_x,area_y,area_p,...
           nsup_x,nsup_y,nsup_p,hx_x,hx_y,hx_p,hy_x,hy_y,hy_p,theta);
     end
     if imovie==1
       dummy=movieout(icount,xc,yc,rc,xp,yp,xlag,ylag,u,v);
     end
  end
% update the lagrangian nodes and related variables
%   TORQUE=sum((xlag(nlsts:nlste)-xlag(nlsts)).*FY'-(ylag(nlsts:nlste)-ylag(nlsts)).*FX');
   TORQUE=sum((xlag(nlsts:nlste)-xc).*FY'-(ylag(nlsts:nlste)-yc).*FX');
   if mod(it,iconv) == 0
     disp('Torque'); disp(TORQUE);
     plot(xlag,ylag,'k',xc,yc,'r*');
     drawnow;
   end
if iflap==1
   rhs_nm1=coef_rot1*angolo-coef_rot2*TORQUE;
   if mod(it,iconv) == 0
     disp('Elastic reaction'); disp(-coef_rot1*angolo);
   end
   if AB2==0
     rhs_nm2=0;
   end
   angolo_dot_old=angolo_dot;
   angolo_dot=angolo_dot_old-Deltat*AB1*rhs_nm1+Deltat*AB2*rhs_nm2;
   angolo=angolo+Deltat/2*(angolo_dot_old+angolo_dot);
   rhs_nm2=rhs_nm1;
   [xlag,ylag,epsilonx,epsilony,bcoefx,bcoefy,ix_x,ix_y,jy_x,jy_y,area_x,area_y,nsup_x,nsup_y,hx_x,hx_y,hy_x,hy_y]=...
   update_lag_tail(xu,yu,xv,yv,dS,nlsts,nlste,xlag,ylag,angolo);
%   [xlag,ylag,epsilonx,epsilony,bcoefx,bcoefy,ix_x,ix_y,jy_x,jy_y,area_x,area_y,nsup_x,nsup_y,hx_x,hx_y,hy_x,hy_y]=...
%   update_lag_onlytail(xu,yu,xv,yv,dS,nlsts,nlste,xlag,ylag,xc,yc,center_rot,angolo);
    if mod(it,iconv) == 0
     disp('angolo'); disp(angolo);
     disp('angolo_dot'); disp(angolo_dot);
   end
end
toc
disp('time for 1 time step')
AB1=3/2; AB2=1/2; CGP=1; % Adam Bashfort coef.
end % END of time looping


%%%%%%%%%%%%%%%%%%%%%% ENDING STORY


%%% some post process
[omega,psi,xo,yo]=vort(u,v,hx,hy);
[X,Y]=meshgrid(xo(2:ncx),yo(2:ncy));
%[omega,psi,xo,yo]=vort(u,v,hx,hy);
fp=figure;
ap=axes;
maxval=max(max(psi));
minval=min(min(psi));
vals=linspace(minval,maxval,50);
contour(X,Y,psi',vals);
set(ap,'DataAspectRatio',[1 1 1],'XLim',[0 lx],'YLim',[0 ly]);
title(sprintf(' Stream function  (%g...%g); Re_b=%g',...
                minval,maxval,Re));
colorbar;

fo=figure;
ao=axes;
maxval=max(max(omega));
minval=min(min(omega));
vals=linspace(minval,maxval,100);
contour(X,Y,omega(2:ncx,2:ncy)',vals);
set(ao,'DataAspectRatio',[1 1 1],'XLim',[0 lx],'YLim',[0 ly]);
title(sprintf(' Vorticity (%g...%g); Re_b=%g',...
                minval,maxval,Re));
colorbar;

fx=figure;
ax=axes;
maxval=max(max(u));
minval=min(min(u));
vals=linspace(minval,maxval,50);
contour(xu,yu,u',vals);
set(ax,'DataAspectRatio',[1 1 1],'XLim',[0 lx],'YLim',[0 ly]);
title(sprintf(' u-com velocity (%g...%g); Re_b=%g',...
                minval,maxval,Re));
colorbar;


