clear
close all

makeplot=1;
makevid=0;
overlay=0;

m=10; %mass
D=1; %drag coefficient
g=9.8;
omega=1;
A0=1;
gamma=A0*omega^2;
k=1; %spring
b=1; %damping

dt=0.1; %time step
t0=0; %initial time
Nt=300; %number of timesteps after initial time t0 (total = Nt+1)
Ntt=100; %number of interior points to evaluate in ode45
tmax=(Nt)*dt+t0; %end time
tvals=t0:dt:tmax; %time values
ttvals=t0:dt/Ntt:tmax; %time values including interior points

dx=0.1; %space step
x0=0; %center
Nx=600; %number of space steps left and right (total = 2Nx+1)
xmin=x0-Nx*dx; %leftmost
xmax=x0+Nx*dx; %rightmost
xvals=xmin:dx:xmax; %all space values

bouncedt=1; %time between bounces
droptimes=t0:bouncedt:tmax; %times when drop occurs
numdrops=floor(tmax/bouncedt)+1; %number of times droplet hits surface
dropposns=zeros(1,numdrops); %positions where drop occurs
dropcount=0; %drops so far

H = zeros(length(tvals),length(xvals)); % wave field
gradH = zeros(length(tvals),length(xvals)); % gradient of wave field
full_trajectory = zeros(1,(Nt)*Ntt+1);
full_velocities = zeros(1,(Nt)*Ntt+1);

% wave created by one drop:
h = @(x0,x,t) besselj(0,omega*abs(x-x0)) * exp(-t); %x0=point of contact


if makevid==1
    vidObj = VideoWriter('WaveVideo7.mp4','MPEG-4');
    vidObj.FrameRate = 10;
    vidObj.Quality = 100;
    open(vidObj)
end
if makeplot==1 || makevid==1
    wave=figure(1);
    clf
end

%initialize x,xprime
x=x0;
xprime=0;

% loop through iterations of time
for i=1:Nt+1
    t=tvals(i);

    %check if drop occurs
    if dropcount<numdrops && t==droptimes(dropcount+1)
        dropposns(dropcount+1)=x;
        dropcount=dropcount+1;
    end

    % calculate wave field
    for j=1:dropcount
        % wave field at time tvals(i)=t
        H(i,:) = H(i,:) + h(dropposns(j),xvals,t-droptimes(j));
    end

    % gradient of wave field
    gradH(i,:) = gradient(H(i,:));

    % solve for trajectory of drop between this timestep (t) and next
    % timestep (t+dt)
    if i<Nt+1
        r=2*rand(1,2)-1;
        x=x+r(1)*0.01;
        xprime=xprime+r(2)*0.01;
        x_ind = find(abs(xvals-x)<=dx/2,1);
        tspan=t:dt/Ntt:t+dt;
        x_init=[x; xprime];
        [tnew,xnew] = ode45(@(tt,xx) trajectory_dxdt(tt,xx,gradH(i,x_ind),m,g,D),tspan,x_init);
        % tnew: Ntt+1 by 1; xnew: Ntt+1 by 2

        full_trajectory((i-1)*Ntt+1:i*Ntt) = xnew(1:Ntt,1);
        full_velocities((i-1)*Ntt+1:i*Ntt) = xnew(1:Ntt,2);
    else
        full_trajectory(end) = xnew(end,1);
        full_velocities(end) = xnew(end,2);
    end

    %set new x and xprime for next timestep (t+dt)
    x=xnew(end,1);
    xprime=xnew(end,2);

    % plot wave
    if makeplot==1 || makevid==1
        if overlay==1
            hold on
        end
        plot(xvals,0*xvals,'black','Linewidth',1)
        hold on
        plot(xvals,H(i,:),'blue','LineWidth',2)
        plot(xvals,gradH(i,:),'red','LineWidth',2)
        scatter(dropposns(1:dropcount),0*(1:dropcount),'black','filled')
        hold off
        xlabel('horizontal position')
        ylabel('height of wave')
        title('Wave Field')
        xlim([xmin,xmax])
        ylim([-1,2])
        legend('x=0','wave field','gradient','contact point')
        if makevid==1
            frame = getframe(wave);
            writeVideo(vidObj,frame);
        end
    end
end

% plot full trajectory
if makeplot==1
    figure(2)
    clf
    hold on
    plot(ttvals,full_trajectory,'black','LineWidth',2)
    plot(ttvals,full_velocities,'red','LineWidth',2)
    hold off
    xlabel('time')
    title('Droplet Trajectory')
    legend(['position';'velocity'])
    ylabel('horizontal position of droplet')
end

if makevid==1
    close(vidObj)
end


zbar = @(z,H) z-H;

Z = @(zb,zbprime,k,b) tanh(100*(-zb))*max(-k*zb -b*zbprime,0);

function dxdt = trajectory_dxdt(~,xx,gradH,m,g,D)
r=2*rand(1,1)-1; %random number in [-1,1]
dxdt = [xx(2); -g*(gradH+r*0.01) - D*xx(2)*(1/m)];
end

