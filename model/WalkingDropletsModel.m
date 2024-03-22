clear
close all

makeplot=1;
makevid=0;
overlay=0;
savefigures=0;

savelocation='./Plots_and_Videos/TwoDroplets';
setnum=7;

MAX_LENGTH=1000;

m=0.1; %mass
D=0.01; %drag coefficient
g=10;
omega=1;
A0=1;
gamma=A0*omega^2;
k=10; %spring (5)
b=5; %damping (1)

dt=0.1; %time step
t0=0; %initial time
Nt=10; %number of timesteps after initial time t0 (total = Nt+1)
Ntt=100; %number of interior points to evaluate in ode45
tmax=(Nt)*dt+t0; %end time
tvals=t0:dt:tmax; %time values
ttvals=t0:dt/Ntt:tmax; %time values including interior points

dx=0.1; %space step
x0=0; %center
Nx=200; %number of space steps left and right (total = 2Nx+1)
xmin=x0-Nx*dx; %leftmost
xmax=x0+Nx*dx; %rightmost
xvals=xmin:dx:xmax; %all space values

droptimes1=zeros(1,MAX_LENGTH); %times when drop occurs
dropxposns1=zeros(1,MAX_LENGTH); %x positions where drop occurs
dropzposns1=zeros(1,MAX_LENGTH); %z positions where drop occurs
dropcount1=0; %drops so far

droptimes2=zeros(1,MAX_LENGTH); %times when drop occurs
dropxposns2=zeros(1,MAX_LENGTH); %x positions where drop occurs
dropzposns2=zeros(1,MAX_LENGTH); %z positions where drop occurs
dropcount2=0; %drops so far

H = zeros(length(tvals),length(xvals)); % wave field: H(t,x)=height of wave at time t and horizontal position x 
gradH = zeros(length(tvals),length(xvals)); % gradient of wave field
full_xtrajectory1 = zeros(1,(Nt)*Ntt+1);
full_xvelocities1 = zeros(1,(Nt)*Ntt+1);
full_ztrajectory1 = zeros(1,(Nt)*Ntt+1);
full_zvelocities1 = zeros(1,(Nt)*Ntt+1);
full_xtrajectory2 = zeros(1,(Nt)*Ntt+1);
full_xvelocities2 = zeros(1,(Nt)*Ntt+1);
full_ztrajectory2 = zeros(1,(Nt)*Ntt+1);
full_zvelocities2 = zeros(1,(Nt)*Ntt+1);

% wave created by one drop:
h = @(x0,xmin,xmax,x,t) (besselj(0,omega*abs(x-x0)) ...
    + besselj(0,omega*abs(x-x0+(xmax-xmin))) ...
    + besselj(0,omega*abs(x-x0-(xmax-xmin)))) * exp(-t); %x0=point of contact


zbar = @(z,H) z-H;


if makevid==1
    vidObj = VideoWriter(strcat(savelocation,'_Vid',int2str(setnum),'.mp4'),'MPEG-4');
    vidObj.FrameRate = 10;
    vidObj.Quality = 100;
    open(vidObj)
end
if makeplot==1 || makevid==1
    wave=figure(1);
    clf
end


% dropcount=1;
% dropposns(1)=x0;

%initialize x,xprime
x1=-1;
xprime1=0;
x2=2;
xprime2=0;

%initialize z,zprime
z1=0.4;
zprime1=0;
z2=0.6;
zprime2=0;

% loop through iterations of time
for i=1:Nt+1
    t=tvals(i);

    % calculate wave field
    for j=1:dropcount1
        % wave field at time tvals(i)=t
        H(i,:) = H(i,:) + h(dropxposns1(j),xmin,xmax,xvals,t-droptimes1(j));
    end
    for j=1:dropcount2
        H(i,:) = H(i,:) + h(dropxposns2(j),xmin,xmax,xvals,t-droptimes2(j));
    end

    % gradient of wave field
    gradH(i,:) = gradient(H(i,:));

    % time derivative of wave field up until this time
    Hprime = (gradient((H(1:i,:))'))';

    % make sure droplet doesn't pass walls
    if x1>=xmax && xprime1>0
        x1=xmax;
        xprime1=-xprime1;
    elseif x1<=xmin && xprime1<0
        x1=xmin;
        xprime1=-xprime1;
    end
    if x2>=xmax && xprime2>0
        x2=xmax;
        xprime2=-xprime2;
    elseif x2<=xmin && xprime2<0
        x2=xmin;
        xprime2=-xprime2;
    end

    % solve for horizontal trajectory of drop between this timestep (t) and
    % next timestep (t+dt)
    if i<Nt+1
        if dropcount1>=1
            %droplet 1
            % r=2*rand(1,2)-1;
            % x1=x1+r(1)*0.01;
            % xprime1=xprime1+r(2)*0.01;
            lastdroppos1 = dropxposns1(dropcount1);
            x_ind1 = find(abs(xvals-lastdroppos1)<=dx/2+0.001,1);
            tspan1=t:dt/Ntt:t+dt;
            x_init1=[x1; xprime1];
            [tnew1,xnew1] = ode45(@(tt,xx) horizontal_trajectory(tt,xx,gradH(i,x_ind1),m,g,D),tspan1,x_init1);
            % tnew: Ntt+1 by 1; xnew: Ntt+1 by 2

            full_xtrajectory1((i-1)*Ntt+1:i*Ntt) = xnew1(1:Ntt,1);
            full_xvelocities1((i-1)*Ntt+1:i*Ntt) = xnew1(1:Ntt,2);

            %set new x and xprime for next timestep (t+dt)
            nextx1=xnew1(end,1);
            nextxprime1=xnew1(end,2);
            
        else
            if xprime1==0
                full_xtrajectory1((i-1)*Ntt+1:i*Ntt) = x1;
                full_xvelocities1((i-1)*Ntt+1:i*Ntt) = xprime1;
            else
                full_xtrajectory1((i-1)*Ntt+1:i*Ntt) = (x1+xprime1*dt/Ntt):(xprime1*dt/Ntt):(x1+xprime1*dt);
                full_xvelocities1((i-1)*Ntt+1:i*Ntt) = xprime1*ones(Ntt,1);
            end

            nextx1=x1+xprime1*dt;
            nextxprime1=xprime1;
        end
        if dropcount2>=1
            %droplet 2
            % r=2*rand(1,2)-1;
            % x2=x2+r(1)*0.01;
            % xprime2=xprime2+r(2)*0.01;
            lastdroppos2 = dropxposns2(dropcount2);
            x_ind2 = find(abs(xvals-lastdroppos2)<=dx/2+0.001,1);
            tspan2=t:dt/Ntt:t+dt;
            x_init2=[x2; xprime2];
            [tnew2,xnew2] = ode45(@(tt,xx) horizontal_trajectory(tt,xx,gradH(i,x_ind2),m,g,D),tspan2,x_init2);

            full_xtrajectory2((i-1)*Ntt+1:i*Ntt) = xnew2(1:Ntt,1);
            full_xvelocities2((i-1)*Ntt+1:i*Ntt) = xnew2(1:Ntt,2);

            nextx2=xnew2(end,1);
            nextxprime2=xnew2(end,2);
        else
            if xprime2==0
                full_xtrajectory2((i-1)*Ntt+1:i*Ntt) = x2;
                full_xvelocities2((i-1)*Ntt+1:i*Ntt) = xprime2;
            else
                full_xtrajectory2((i-1)*Ntt+1:i*Ntt) = (x2+xprime2*dt/Ntt):(xprime2*dt/Ntt):(x2+xprime2*dt);
                full_xvelocities2((i-1)*Ntt+1:i*Ntt) = xprime2*ones(Ntt,1);
            end

            nextx2=x2+xprime2*dt;
            nextxprime2=xprime2;
        end
    
    
    else
        full_xtrajectory1(end) = x1;
        full_xvelocities1(end) = xprime1;
        full_xtrajectory2(end) = x2;
        full_xvelocities2(end) = xprime2;

        nextx1=x1;
        nextxprime1=xprime1;
        nextx2=x2;
        nextxprime2=xprime2;
    end

    

    % solve for vertical trajectory of drop between this timestep (t) and
    % next timestep (t+dt)
    if i<Nt+1
        x_ind1 = find(abs(xvals-x1)<=dx/2+0.001,1);
        tspan1=t:dt/Ntt:t+dt;
        zb1 = zbar(z1,H(i,x_ind1));
        z_init1 = [z1; zprime1; zb1];
        [tnew1,znew1] = ode45(@(tt,zz) vertical_trajectory(tt,zz,Hprime(i,x_ind1),m,g,gamma,k,b),tspan1,z_init1);

        x_ind2 = find(abs(xvals-x2)<=dx/2+0.001,1);
        tspan2=t:dt/Ntt:t+dt;
        zb2 = zbar(z2,H(i,x_ind2));
        z_init2 = [z2; zprime2; zb2];
        [tnew2,znew2] = ode45(@(tt,zz) vertical_trajectory(tt,zz,Hprime(i,x_ind2),m,g,gamma,k,b),tspan2,z_init2);

        full_ztrajectory1((i-1)*Ntt+1:i*Ntt) = znew1(1:Ntt,1);
        full_zvelocities1((i-1)*Ntt+1:i*Ntt) = znew1(1:Ntt,2);

        full_ztrajectory2((i-1)*Ntt+1:i*Ntt) = znew2(1:Ntt,1);
        full_zvelocities2((i-1)*Ntt+1:i*Ntt) = znew2(1:Ntt,2);

        %set new z and zprime for next timestep (t+dt)
        nextz1=znew1(end,1);
        nextzprime1=znew1(end,2);
        nextz2=znew2(end,1);
        nextzprime2=znew2(end,2);
    else
        full_ztrajectory1(end) = znew1(end,1);
        full_zvelocities1(end) = znew1(end,2);
        full_ztrajectory2(end) = znew2(end,1);
        full_zvelocities2(end) = znew2(end,2);

        nextz1=z1;
        nextzprime1=zprime1;
        nextz2=z2;
        nextzprime2=zprime2;
    end
    

    %check if drop occurs
    if z1>=H(i,x_ind1) && nextz1<=H(i,x_ind1)
        dropxposns1(dropcount1+1)=x1;
        dropzposns1(dropcount1+1)=z1;
        droptimes1(dropcount1+1)=t;
        dropcount1=dropcount1+1;
    end
    if z2>=H(i,x_ind2) && nextz2<=H(i,x_ind2)
        dropxposns2(dropcount2+1)=x2;
        dropzposns2(dropcount2+1)=z2;
        droptimes2(dropcount2+1)=t;
        dropcount2=dropcount2+1;
    end

    %set new positions and velocities
    x1=nextx1;
    xprime1=nextxprime1;
    z1=nextz1;
    zprime1=nextzprime1;
    x2=nextx2;
    xprime2=nextxprime2;
    z2=nextz2;
    zprime2=nextzprime2;
    

    % plot wave
    if makeplot==1 || makevid==1
        if overlay==1
            hold on
        end
        plot(xvals,0*xvals,'black','Linewidth',1)
        hold on
        plot(xvals,H(i,:),'blue','LineWidth',2)
        plot(xvals,gradH(i,:),'red','LineWidth',2)
        scatter([x1,x2],[z1,z2],'blue','filled')
        hold off
        xlabel('horizontal position')
        ylabel('height of wave')
        title('Wave Field')
        xlim([xmin,xmax])
        ylim([-1,2])
        legend('x=0','wave field','gradient','droplets')
        if makevid==1
            frame = getframe(wave);
            writeVideo(vidObj,frame);
        end
    end
end

% plot final wave field
if makeplot==1
    plot(xvals,0*xvals,'black','Linewidth',1)
        hold on
        plot(xvals,H(end,:),'blue','LineWidth',2)
        plot(xvals,gradH(end,:),'red','LineWidth',2)
        scatter([x1,x2],[z1,z2],'blue','filled')
        hold off
        xlabel('horizontal position')
        ylabel('height of wave')
        title('Wave Field')
        xlim([xmin,xmax])
        ylim([-1,2])
        legend('x=0','wave field','gradient','droplets')
end

% plot full horizontal trajectory
if makeplot==1
    horizontal=figure(2);
    clf
    hold on
    plot(ttvals,full_xtrajectory1,'black','LineWidth',2)
    plot(ttvals,full_xvelocities1,'red','LineWidth',2)
    plot(ttvals,full_xtrajectory2,'black','LineWidth',2,'LineStyle','--')
    plot(ttvals,full_xvelocities2,'red','LineWidth',2,'LineStyle','--')
    hold off
    xlabel('time')
    title('Horizontal Trajectory')
    legend('droplet 1 position','droplet 1 velocity','droplet 2 position','droplet 2 velocity')
    ylabel('horizontal position of droplet')
    if savefigures==1
        savefig(horizontal,strcat(savelocation,'_Horizontal',int2str(setnum),'.fig'))
    end
end

% plot full vertical trajectory
if makeplot==1
    vertical=figure(3);
    clf
    hold on
    plot(ttvals,0*ttvals,'black','Linewidth',1)
    plot(ttvals,full_ztrajectory1,'black','LineWidth',2)
    plot(ttvals,full_zvelocities1,'red','LineWidth',2)
    plot(ttvals,full_ztrajectory2,'black','LineWidth',2,'LineStyle', '--')
    plot(ttvals,full_zvelocities2,'red','LineWidth',2,'LineStyle', '--')
    hold off
    xlabel('time')
    title('Vertical Trajectory')
    legend('x=0','droplet 1 position','droplet 1 velocity','droplet 2 position','droplet 2 velocity')
    ylabel('vertical position of droplet')
    if savefigures==1
        savefig(vertical,strcat(savelocation,'_Vertical',int2str(setnum),'.fig'))
    end
end

if makevid==1
    close(vidObj)
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

