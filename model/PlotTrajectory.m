close all
clear

load('./Plots_and_Videos/2D/MatchingInitialPositions_0_5Mass8_data.mat')
% load("Plots_and_Videos/2D/TwoDroplets11_data.mat")

timeskip=250;
speedthreshold=0.5;

red=[1,0,0];
blue=[0,0,1];
colorLength=100;
redtoblue=[linspace(red(1),blue(1),colorLength)',linspace(red(2),blue(2),colorLength)',linspace(red(3),blue(3),colorLength)'];

m=0.5; %mass
D=0.01; %drag coefficient
g=10;
omega=1;
A0=1;
gamma=A0*omega^2;
k=10; %spring (5)
b=0.1; %damping (1)

dt=0.1; %time step
t0=0; %initial time
% Nt=2000; %number of timesteps after initial time t0 (total = Nt+1)
Nt=(length(full_x1trajectory)-1)/100;
Ntt=100; %number of interior points to evaluate in ode45
tmax=(Nt)*dt+t0; %end time
tvals=t0:dt:tmax; %time values
ttvals=t0:dt/Ntt:tmax; %time values including interior points

dx=0.1; %x step
x0=0; %center
Nx=150; %number of space steps left and right (total = 2Nx+1)
xmin=x0-Nx*dx; %leftmost
xmax=x0+Nx*dx; %rightmost
xvals=xmin:dx:xmax; %all x values

dy=0.1; %y step
y0=0; %center
Ny=150; %number of y steps up and down (total = 2Ny+1)
ymin=y0-Ny*dy; %lowest
ymax=y0+Ny*dy; %highest
yvals=ymin:dy:ymax; %all y values

boundary_radius=xmax-x0;

horizontal=figure(1);
clf
hold on
speedvals1=sqrt(full_x1trajectory(2,:).^2+full_y1trajectory(2,:).^2);
if max(speedvals1>0)
    speedvals1=speedvals1/max(speedvals1);
end
speedvals2=sqrt(full_x2trajectory(2,:).^2+full_y2trajectory(2,:).^2);
if max(speedvals2>0)
    speedvals2=speedvals2/max(speedvals2);
end
speedvals3=sqrt(full_x3trajectory(2,:).^2+full_y3trajectory(2,:).^2);
if max(speedvals3>0)
    speedvals3=speedvals3/max(speedvals3);
end
speedvals4=sqrt(full_x4trajectory(2,:).^2+full_y4trajectory(2,:).^2);
if max(speedvals4>0)
    speedvals4=speedvals4/max(speedvals4);
end
speedvals5=sqrt(full_x5trajectory(2,:).^2+full_y5trajectory(2,:).^2);
if max(speedvals5>0)
    speedvals5=speedvals5/max(speedvals5);
end
scolors1=speedvals1'*red+(1-speedvals1')*blue;
scolors2=speedvals2'*red+(1-speedvals2')*blue;
scolors3=speedvals3'*red+(1-speedvals3')*blue;
scolors4=speedvals4'*red+(1-speedvals4')*blue;
scolors5=speedvals5'*red+(1-speedvals5')*blue;
alphas1=speedvals1;
alphas1(alphas1<speedthreshold)=alphas1(alphas1<speedthreshold)/5;
alphas2=speedvals2;
alphas2(alphas2<speedthreshold)=alphas2(alphas2<speedthreshold)/5;
alphas3=speedvals3;
alphas3(alphas3<speedthreshold)=alphas3(alphas3<speedthreshold)/5;
alphas4=speedvals4;
alphas4(alphas4<speedthreshold)=alphas4(alphas4<speedthreshold)/5;
alphas5=speedvals5;
alphas5(alphas5<speedthreshold)=alphas5(alphas5<speedthreshold)/5;
scatter(full_x1trajectory(1,1:timeskip:end),full_y1trajectory(1,1:timeskip:end),[],scolors1(1:timeskip:end,:),'filled','AlphaData',alphas1(1:timeskip:end),'MarkerFaceAlpha','flat');
scatter(full_x2trajectory(1,1:timeskip:end),full_y2trajectory(1,1:timeskip:end),[],scolors2(1:timeskip:end,:),'filled','AlphaData',alphas2(1:timeskip:end),'MarkerFaceAlpha','flat');
scatter(full_x3trajectory(1,1:timeskip:end),full_y3trajectory(1,1:timeskip:end),[],scolors3(1:timeskip:end,:),'filled','AlphaData',alphas3(1:timeskip:end),'MarkerFaceAlpha','flat');
scatter(full_x4trajectory(1,1:timeskip:end),full_y4trajectory(1,1:timeskip:end),[],scolors4(1:timeskip:end,:),'filled','AlphaData',alphas4(1:timeskip:end),'MarkerFaceAlpha','flat');
scatter(full_x5trajectory(1,1:timeskip:end),full_y5trajectory(1,1:timeskip:end),[],scolors5(1:timeskip:end,:),'filled','AlphaData',alphas5(1:timeskip:end),'MarkerFaceAlpha','flat');
% surface([full_x1trajectory(1,:);full_x1trajectory(1,:)],[full_y1trajectory(1,:);full_y1trajectory(1,:)],...
%     [zeros(size(full_x1trajectory(1,:)));zeros(size(full_x1trajectory(1,:)))],[speedvals1;speedvals1],...
%     'facecol','no',...
%     'edgecol','interp',...
%     'linew',2);
% surface([full_x2trajectory(1,:);full_x2trajectory(1,:)],[full_y2trajectory(1,:);full_y2trajectory(1,:)],...
%     [zeros(size(full_x2trajectory(1,:)));zeros(size(full_x2trajectory(1,:)))],[speedvals2;speedvals2],...
%     'facecol','no',...
%     'edgecol','interp',...
%     'linew',2);
% s.MarkerFaceAlpha='flat';
% colormap('cool')
colormap(flip(redtoblue))
cb=colorbar;
% cb.Label.String = 'time';
% clim([ttvals(1),ttvals(end)])
cb.Label.String = 'Normalized Speed';
clim([0,1])
% xlim([xmin, xmax])
% ylim([ymin, ymax])
scatter3([full_x1trajectory(1,1),full_x2trajectory(1,1),full_x3trajectory(1,1),full_x4trajectory(1,1),full_x5trajectory(1,1)],...
    [full_y1trajectory(1,1),full_y2trajectory(1,1),full_y3trajectory(1,1),full_y4trajectory(1,1),full_y5trajectory(1,1)],...
    zeros(1,5),100,'black','filled')
scatter3([full_x1trajectory(1,end),full_x2trajectory(1,end),full_x3trajectory(1,end),full_x4trajectory(1,end),full_x5trajectory(1,end)],...
    [full_y1trajectory(1,end),full_y2trajectory(1,end),full_y3trajectory(1,end),full_y4trajectory(1,end),full_y5trajectory(1,end)],...
    zeros(1,5),100,'black','filled','square')
% scatter3([full_x1trajectory(1,1),full_x2trajectory(1,1)],...
%     [full_y1trajectory(1,1),full_y2trajectory(1,1)],...
%     zeros(1,2),100,'black','filled')
% scatter3([full_x1trajectory(1,end),full_x2trajectory(1,end)],...
%     [full_y1trajectory(1,end),full_y2trajectory(1,end)],...
%     zeros(1,2),100,'black','filled','square')
% plot(full_xtrajectory,full_ytrajectory,'MarkerFaceColor',(ttvals/max(ttvals))'*[1,0,0]+(1-ttvals/max(ttvals))'*[0,0,1],'LineWidth',2)
% plot(ttvals,full_xvelocities,'red','LineWidth',2)
th=0:pi/50:2*pi;
plot(boundary_radius*cos(th),boundary_radius*sin(th),'black')
hold off
xlabel('x')
title('Modeled Trajectories')
% legend(['position';'start';'end'])
legend('','','','','','start','end')
ylabel('y')
xlim([xmin,xmax])
ylim([ymin,ymax])
set(gcf,'color','white')





% figure(2)
% clf
% hold on
% plot(ttvals,0*ttvals,'black','Linewidth',1)
% plot(ttvals,full_z1trajectory(1,:),'black','LineWidth',2)
% plot(ttvals,full_z1trajectory(2,:),'red','LineWidth',2)
% plot(ttvals,full_z2trajectory(1,:),'black','LineWidth',2,'LineStyle','--')
% plot(ttvals,full_z2trajectory(2,:),'red','LineWidth',2,'LineStyle','--')
% hold off
% xlabel('time')
% title('Vertical Trajectory')
% legend('x=0','droplet 1 position','droplet 1 velocity','droplet 2 position', 'droplet 2 velocity')
% ylabel('vertical position of droplet')
% set(gcf,'color','white')