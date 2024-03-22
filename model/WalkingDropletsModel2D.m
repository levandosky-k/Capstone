clear
close all

makeplot=1;
makevid=1;
overlay=0;
savefigures=1;

savelocation='./Plots_and_Videos/2D/InitialXYVelocity';
setnum=1;

MAX_LENGTH=1000;

m=0.2; %mass
D=0.01; %drag coefficient
g=10;
omega=1;
A0=1;
gamma=A0*omega^2;
k=10; %spring (5)
b=0.1; %damping (1)

dt=0.1; %time step
t0=0; %initial time
Nt=300; %number of timesteps after initial time t0 (total = Nt+1)
Ntt=100; %number of interior points to evaluate in ode45
tmax=(Nt)*dt+t0; %end time
tvals=t0:dt:tmax; %time values
ttvals=t0:dt/Ntt:tmax; %time values including interior points

dx=0.1; %x step
x0=0; %center
Nx=100; %number of space steps left and right (total = 2Nx+1)
xmin=x0-Nx*dx; %leftmost
xmax=x0+Nx*dx; %rightmost
xvals=xmin:dx:xmax; %all x values

dy=0.1; %y step
y0=0; %center
Ny=100; %number of y steps up and down (total = 2Ny+1)
ymin=y0-Ny*dy; %lowest
ymax=y0+Ny*dy; %highest
yvals=ymin:dy:ymax; %all y values

bouncedt=10; %time between bounces
numdrops=floor(tmax/bouncedt)+1; %number of times droplet hits surface
% droptimes=t0:bouncedt:tmax; %times when drop occurs
droptimes=zeros(1,MAX_LENGTH);
dropxposns=zeros(1,MAX_LENGTH); %x positions where drop occurs
dropyposns=zeros(1,MAX_LENGTH); %y positions where drop occurs
dropzposns=zeros(1,MAX_LENGTH); %z positions where drop occurs
dropcount=0; %drops so far

H = zeros(length(xvals)*length(tvals),length(yvals)); % wave field
gradH_X = zeros(length(xvals)*length(tvals),length(yvals)); % gradient of wave field in x direction
gradH_Y = zeros(length(xvals)*length(tvals),length(yvals));
full_xtrajectory = zeros(1,(Nt)*Ntt+1);
full_xvelocities = zeros(1,(Nt)*Ntt+1);
full_ytrajectory = zeros(1,(Nt)*Ntt+1);
full_yvelocities = zeros(1,(Nt)*Ntt+1);
full_ztrajectory = zeros(1,(Nt)*Ntt+1);
full_zvelocities = zeros(1,(Nt)*Ntt+1);

% wave created by one drop:
% h = @(x0,y0,x,y,t) besselj(0,omega*abs([x;y]-[x0;y0])) * exp(-t); %(x0,y0)=point of contact

zbar = @(z,H) z-H;


if makevid==1
    vidObj = VideoWriter(strcat(savelocation,'_Vid',int2str(setnum),'.mp4'),'MPEG-4');
    vidObj.FrameRate = 20;
    vidObj.Quality = 100;
    open(vidObj)
end
if makeplot==1 || makevid==1
    wave=figure(1);
    clf
end

% dropcount=1;
% dropposns(1)=x0;

%initialize x,xprime,y,yprime
x=0;
xprime=1;
y=0;
yprime=2;

%initialize z,zprime
z=0.7;
zprime=0;

hit=0;

tic
% loop through iterations of time
for i=1:Nt+1
    t=tvals(i);

    %check if drop occurs
    % if dropcount<numdrops && t==droptimes(dropcount+1)
    %     dropxposns(dropcount+1)=x;
    %     dropyposns(dropcount+1)=y;
    %     dropcount=dropcount+1;
    % end

    % calculate wave field
    for j=1:dropcount
        % wave field at time tvals(i)=t
        H(((i-1)*length(xvals)+1):(i*length(xvals)),:) = H(((i-1)*length(xvals)+1):(i*length(xvals)),:) + single_wave(dropxposns(j),dropyposns(j),xvals,yvals,t-droptimes(j),omega);
    end

    % size(H(i*(1:length(xvals)),:)==0)
    % H(((i-1)*length(xvals)+1):(i*length(xvals)),:) = H(((i-1)*length(xvals)+1):(i*length(xvals)),:) + single_wave(dropxposns(1),dropyposns(1),xvals,yvals,t-droptimes(1),omega);

    % gradient of wave field at this timestep
    [FY,FX] = gradient(H(((i-1)*length(xvals)+1):(i*length(xvals)),:));
    gradH_X(((i-1)*length(xvals)+1):(i*length(xvals)),:) = FX;
    gradH_Y(((i-1)*length(xvals)+1):(i*length(xvals)),:) = FY;
    
    

    % make sure droplet doesn't pass walls
    if x>=xmax && xprime>0
        x=xmax;
        xprime=-xprime;
    elseif x<=xmin && xprime<0
        x=xmin;
        xprime=-xprime;
    end
    if y>=ymax && yprime>0
        y=ymax;
        yprime=-yprime;
    elseif y<=ymin && yprime<0
        y=ymin;
        yprime=-yprime;
    end


    % solve for horizontal trajectory of drop between this timestep (t) and
    % next timestep (t+dt)
    if i<Nt+1
        if dropcount>=1
            lastdropposx = dropxposns(dropcount);
            lastdropposy = dropyposns(dropcount);
            x_ind = find(abs(xvals-lastdropposx)<=dx/2+0.001,1);
            y_ind = find(abs(yvals-lastdropposy)<=dx/2+0.001,1);
            % x_ind = find(abs(xvals-x)<=dx/2,1);
            % y_ind = find(abs(yvals-y)<=dy/2,1);
            tspan=t:dt/Ntt:t+dt;
            x_init=[x; xprime];
            y_init=[y; yprime];
            [tnew,xnew] = ode45(@(tt,xx) horizontal_trajectory(tt,xx,gradH_X(i*x_ind,y_ind),m,g,D),tspan,x_init);
            [~,ynew] = ode45(@(tt,yy) horizontal_trajectory(tt,yy,gradH_Y(i*x_ind,y_ind),m,g,D),tspan,y_init);
            % tnew: Ntt+1 by 1; xnew: Ntt+1 by 2

            full_xtrajectory((i-1)*Ntt+1:i*Ntt) = xnew(1:Ntt,1);
            full_xvelocities((i-1)*Ntt+1:i*Ntt) = xnew(1:Ntt,2);
            full_ytrajectory((i-1)*Ntt+1:i*Ntt) = ynew(1:Ntt,1);
            full_yvelocities((i-1)*Ntt+1:i*Ntt) = ynew(1:Ntt,2);

            nextx=xnew(end,1);
            nextxprime=xnew(end,2);
            nexty=ynew(end,1);
            nextyprime=ynew(end,2);
        else
            if xprime==0
                full_xtrajectory((i-1)*Ntt+1:i*Ntt) = x;
                full_xvelocities((i-1)*Ntt+1:i*Ntt) = xprime;
            else
                full_xtrajectory((i-1)*Ntt+1:i*Ntt) = (x+xprime*dt/Ntt):(xprime*dt/Ntt):(x+xprime*dt);
                full_xvelocities((i-1)*Ntt+1:i*Ntt) = xprime*ones(Ntt,1);
            end

            if yprime==0
                full_ytrajectory((i-1)*Ntt+1:i*Ntt) = y;
                full_yvelocities((i-1)*Ntt+1:i*Ntt) = yprime;
            else
                full_ytrajectory((i-1)*Ntt+1:i*Ntt) = (y+yprime*dt/Ntt):(yprime*dt/Ntt):(y+yprime*dt);
                full_yvelocities((i-1)*Ntt+1:i*Ntt) = yprime*ones(Ntt,1);
            end

            nextx=x+xprime*dt;
            nextxprime=xprime;
            nexty=y+yprime*dt;
            nextyprime=yprime;
        end
    else
        full_xtrajectory(end) = x;
        full_xvelocities(end) = xprime;
        full_ytrajectory(end) = y;
        full_yvelocities(end) = yprime;

        nextx=x;
        nextxprime=xprime;
        nexty=y;
        nextyprime=yprime;
    end


    % solve for vertical trajectory of drop between this timestep (t) and
    % next timestep (t+dt)
    if i<Nt+1
        x_ind = find(abs(xvals-x)<=dx/2,1);
        y_ind = find(abs(yvals-y)<=dy/2,1);
        tspan=t:dt/Ntt:t+dt;
        zb = zbar(z,H(i*x_ind,y_ind));
        z_init = [z; zprime; zb];
        % calculate time derivative of wave field up to this time
        Hprime = gradient(H(x_ind:length(xvals):i*length(xvals),y_ind));
        
        [tnew,znew] = ode45(@(tt,zz) vertical_trajectory(tt,zz,Hprime(i),m,g,gamma,k,b),tspan,z_init);

        full_ztrajectory((i-1)*Ntt+1:i*Ntt) = znew(1:Ntt,1);
        full_zvelocities((i-1)*Ntt+1:i*Ntt) = znew(1:Ntt,2);

        nextz=znew(end,1);
        nextzprime=znew(end,2);
    else
        full_ztrajectory(end) = znew(end,1);
        full_zvelocities(end) = znew(end,2);

        nextz=z;
        nextzprime=zprime;
    end

    
    x_ind = find(abs(xvals-x)<=dx/2,1);
    y_ind = find(abs(yvals-y)<=dy/2,1);

    % check if drop occurs
    if z>=max(H(i*x_ind,y_ind),0) && nextz<=max(H(i*x_ind,y_ind),0)
    % if z>=0 && nextz<=0
        % hit=1;
        dropxposns(dropcount+1)=x;
        dropyposns(dropcount+1)=y;
        dropzposns(dropcount+1)=z;
        droptimes(dropcount+1)=t;
        dropcount=dropcount+1;
    end

    % if hit==1 && z>=0.2
    %     dropxposns(dropcount+1)=x;
    %     dropzposns(dropcount+1)=z;
    %     droptimes(dropcount+1)=t;
    %     dropcount=dropcount+1;
    %     hit=0;
    % end


    %set new values for next timestep (t+dt)
    x=nextx;
    xprime=nextxprime;
    y=nexty;
    yprime=nextyprime;
    z=nextz;
    zprime=nextzprime;

    % calculate wave field
    % for j=1:dropcount
    %     % wave field at time tvals(i)=t
    %     H(i,:) = H(i,:) + h(dropxposns(j),xvals,t-droptimes(j));
    % end

    % plot wave
    if makevid==1
        if overlay==1
            hold on
        end
        figure(1)
        clf
        [xplot,yplot] = meshgrid(xvals,yvals);
        surf(xplot,yplot,H(((i-1)*length(xvals)+1):((i)*length(xvals)),:))
        alpha 0.5
        % view(2)
        hold on
        scatter3(x,y,z,'red','filled')
        hold on
        shading interp
        colormap('winter')
        cb=colorbar;
        cb.Label.String = 'wave height';
        clim([0,1])
        xlim([xmin,xmax])
        ylim([ymin,ymax])
        zlim([-2,2])
        hold off
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Wave Field')
        if makevid==1
            frame = getframe(wave);
            writeVideo(vidObj,frame);
        end
    end
end
toc

% plot final wave field
% if makeplot==1
%     figure(1)
%     clf
%     [xplot,yplot] = meshgrid(xvals,yvals);
%     surf(xplot,yplot,H(((Nt)*length(xvals)+1):((Nt+1)*length(xvals)),:))
%     view(2)
%     hold on
%     shading interp
%     colormap('winter')
%     cb=colorbar;
%     cb.Label.String = 'wave height';
%     clim([0,1])
%     hold off
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     title('Wave Field')
% end

% plot full x trajectory
if makeplot==1
    xfig=figure(2);
    clf
    hold on
    plot(ttvals,full_xtrajectory,'black','LineWidth',2)
    plot(ttvals,full_xvelocities,'red','LineWidth',2)
    hold off
    xlabel('time')
    title('X Trajectory')
    legend(['position';'velocity'])
    ylabel('x position of droplet')
    if savefigures==1
        savefig(xfig,strcat(savelocation,'_X',int2str(setnum),'.fig'))
    end
end
% plot full y trajectory
if makeplot==1
    yfig=figure(3);
    clf
    hold on
    plot(ttvals,full_ytrajectory,'black','LineWidth',2)
    plot(ttvals,full_yvelocities,'red','LineWidth',2)
    hold off
    xlabel('time')
    title('Y Trajectory')
    legend(['position';'velocity'])
    ylabel('y position of droplet')
    if savefigures==1
        savefig(yfig,strcat(savelocation,'_Y',int2str(setnum),'.fig'))
    end
end
% plot full horizontal trajectory with speed as color
if makeplot==1
    horizontal=figure(4);
    clf
    hold on
    % surface([full_xtrajectory;full_xtrajectory],[full_ytrajectory;full_ytrajectory],...
    %     [zeros(size(full_xtrajectory));zeros(size(full_xtrajectory))],[ttvals;ttvals],...
    %     'facecol','no',...
    %     'edgecol','interp',...
    %     'linew',2);
    speedvals=sqrt(full_xvelocities.^2+full_yvelocities.^2);
    surface([full_xtrajectory;full_xtrajectory],[full_ytrajectory;full_ytrajectory],...
        [zeros(size(full_xtrajectory));zeros(size(full_xtrajectory))],[speedvals;speedvals],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    colormap('cool')
    cb=colorbar;
    % cb.Label.String = 'time';
    % clim([ttvals(1),ttvals(end)])
    cb.Label.String = 'speed';
    clim([min(speedvals),max(speedvals)])
    % xlim([xmin, xmax])
    % ylim([ymin, ymax])
    scatter3(full_xtrajectory(1,1),full_ytrajectory(1,1),0,100,'black','filled')
    scatter3(full_xtrajectory(1,end),full_ytrajectory(1,end),0,100,'black','filled','square')
    % plot(full_xtrajectory,full_ytrajectory,'MarkerFaceColor',(ttvals/max(ttvals))'*[1,0,0]+(1-ttvals/max(ttvals))'*[0,0,1],'LineWidth',2)
    % plot(ttvals,full_xvelocities,'red','LineWidth',2)
    hold off
    xlabel('x')
    title('Horizontal Trajectory')
    % legend(['position';'start';'end'])
    legend('','start','end')
    ylabel('y')
    if savefigures==1
        savefig(horizontal,strcat(savelocation,'_Horizontal',int2str(setnum),'.fig'))
    end
end

% plot full vertical trajectory
if makeplot==1
    vertical=figure(5);
    clf
    hold on
    plot(ttvals,0*ttvals,'black','Linewidth',1)
    plot(ttvals,full_ztrajectory,'black','LineWidth',2)
    plot(ttvals,full_zvelocities,'red','LineWidth',2)
    hold off
    xlabel('time')
    title('Vertical Trajectory')
    legend('x=0','vertical position','vertical velocity')
    ylabel('vertical position of droplet')
    if savefigures==1
        savefig(vertical,strcat(savelocation,'_Vertical',int2str(setnum),'.fig'))
    end
end

if makevid==1
    close(vidObj)
end


function h = single_wave(x0,y0,x,y,t,omega)
h=zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        h(j,i) = besselj(0,omega*sqrt((x(i)-x0)^2 + (y(j)-y0)^2)) * exp(-t);
    end
end
end

function dxdt = horizontal_trajectory(~,xx,gradH,m,g,D)
r=2*rand(1,1)-1; %random number in [-1,1]
dxdt = [xx(2); -g*(gradH) - D*xx(2)*(1/m)];
end


function dzdt = vertical_trajectory(~,zz,Hprime,m,g,gamma,k,b)
Z = @(zb,zbprime,k,b) tanh(100*(-zb))*max(-k*zb -b*zbprime,0);
% Z = @(zb,zbprime,k,b) (1*(-zb>0)+0.5*(-zb==0)+0*(-zb<0))*max(-k*zb -b*zbprime,0);
ZZ=Z(zz(3),zz(2)-Hprime,k,b);
dzdt = [zz(2); ZZ - m*(g+gamma); zz(2)-Hprime];
end
