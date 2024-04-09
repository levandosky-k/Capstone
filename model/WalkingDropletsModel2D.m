clear
close all

makeplot=1;
makevid=1;
overlay=0;
savefigures=0;
savematfile=1;

savelocation='./Plots_and_Videos/2D/MatchingInitialPositions_0_5Mass';
setnum=8;

% Initial positions:
  % 268.0615  111.3026
  % 166.1495  119.6493
  % 152.6409  127.8061
  % 287.9382  107.2358
  % 325.5790  163.4047

MAX_LENGTH=1000;

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
Nt=5000; %number of timesteps after initial time t0 (total = Nt+1)
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

bouncedt=10; %time between bounces
% numdrops=floor(tmax/bouncedt)+1; %number of times droplet hits surface
% droptimes=t0:bouncedt:tmax; %times when drop occurs
droptimes1=zeros(1,MAX_LENGTH);
dropxposns1=zeros(1,MAX_LENGTH); %x positions where drop occurs
dropyposns1=zeros(1,MAX_LENGTH); %y positions where drop occurs
dropzposns1=zeros(1,MAX_LENGTH); %z positions where drop occurs
dropcount1=0; %drops so far
droptimes2=zeros(1,MAX_LENGTH);
dropxposns2=zeros(1,MAX_LENGTH); %x positions where drop occurs
dropyposns2=zeros(1,MAX_LENGTH); %y positions where drop occurs
dropzposns2=zeros(1,MAX_LENGTH); %z positions where drop occurs
dropcount2=0; %drops so far
droptimes3=zeros(1,MAX_LENGTH);
dropxposns3=zeros(1,MAX_LENGTH); %x positions where drop occurs
dropyposns3=zeros(1,MAX_LENGTH); %y positions where drop occurs
dropzposns3=zeros(1,MAX_LENGTH); %z positions where drop occurs
dropcount3=0; %drops so far
droptimes4=zeros(1,MAX_LENGTH);
dropxposns4=zeros(1,MAX_LENGTH); %x positions where drop occurs
dropyposns4=zeros(1,MAX_LENGTH); %y positions where drop occurs
dropzposns4=zeros(1,MAX_LENGTH); %z positions where drop occurs
dropcount4=0; %drops so far
droptimes5=zeros(1,MAX_LENGTH);
dropxposns5=zeros(1,MAX_LENGTH); %x positions where drop occurs
dropyposns5=zeros(1,MAX_LENGTH); %y positions where drop occurs
dropzposns5=zeros(1,MAX_LENGTH); %z positions where drop occurs
dropcount5=0; %drops so far

H = zeros(length(xvals)*length(tvals),length(yvals)); % wave field
gradH_X = zeros(length(xvals)*length(tvals),length(yvals)); % gradient of wave field in x direction
gradH_Y = zeros(length(xvals)*length(tvals),length(yvals));
full_x1trajectory = zeros(2,(Nt)*Ntt+1); %position in first row, velocity in second row
% full_x1velocities = zeros(1,(Nt)*Ntt+1);
full_y1trajectory = zeros(2,(Nt)*Ntt+1);
% full_y1velocities = zeros(1,(Nt)*Ntt+1);
full_z1trajectory = zeros(2,(Nt)*Ntt+1);
% full_z1velocities = zeros(1,(Nt)*Ntt+1);
full_x2trajectory = zeros(2,(Nt)*Ntt+1);
full_y2trajectory = zeros(2,(Nt)*Ntt+1);
full_z2trajectory = zeros(2,(Nt)*Ntt+1);
full_x3trajectory = zeros(2,(Nt)*Ntt+1);
full_y3trajectory = zeros(2,(Nt)*Ntt+1);
full_z3trajectory = zeros(2,(Nt)*Ntt+1);
full_x4trajectory = zeros(2,(Nt)*Ntt+1);
full_y4trajectory = zeros(2,(Nt)*Ntt+1);
full_z4trajectory = zeros(2,(Nt)*Ntt+1);
full_x5trajectory = zeros(2,(Nt)*Ntt+1);
full_y5trajectory = zeros(2,(Nt)*Ntt+1);
full_z5trajectory = zeros(2,(Nt)*Ntt+1);

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
% x1=1;
% xprime1=1;
% y1=2;
% yprime1=-3;
% x2=-1;
% xprime2=0;
% y2=5;
% yprime2=2;
x1=(268.0615-200)*0.1;
xprime1=0;
y1=(111.3026-200)*0.1;
yprime1=0;
x2=(166.1495-200)*0.1;
xprime2=0;
y2=(119.6493-200)*0.1;
yprime2=0;
x3=(152.6409-200)*0.1;
xprime3=0;
y3=(127.8061-200)*0.1;
yprime3=0;
x4=(287.9382-200)*0.1;
xprime4=0;
y4=(107.2358-200)*0.1;
yprime4=0;
x5=(325.5790-200)*0.1;
xprime5=0;
y5=(163.4047-200)*0.1;
yprime5=0;

%initialize z,zprime
z1=0.5;
zprime1=0;
z2=0.5;
zprime2=0;
z3=0.5;
zprime3=0;
z4=0.5;
zprime4=0;
z5=0.5;
zprime5=0;


tic
% loop through iterations of time
for i=1:Nt+1
    t=tvals(i);

    % calculate wave field
    for j=max(1,dropcount1-3):dropcount1
        % wave field at time tvals(i)=t
        H(((i-1)*length(xvals)+1):(i*length(xvals)),:) = H(((i-1)*length(xvals)+1):(i*length(xvals)),:) ...
            + single_wave(dropxposns1(j),dropyposns1(j),xvals,yvals,t-droptimes1(j),omega);
    end
    for j=max(1,dropcount2-3):dropcount2
        % wave field at time tvals(i)=t
        H(((i-1)*length(xvals)+1):(i*length(xvals)),:) = H(((i-1)*length(xvals)+1):(i*length(xvals)),:) ...
            + single_wave(dropxposns2(j),dropyposns2(j),xvals,yvals,t-droptimes2(j),omega);
    end
    for j=max(1,dropcount3-3):dropcount3
        % wave field at time tvals(i)=t
        H(((i-1)*length(xvals)+1):(i*length(xvals)),:) = H(((i-1)*length(xvals)+1):(i*length(xvals)),:) ...
            + single_wave(dropxposns3(j),dropyposns3(j),xvals,yvals,t-droptimes3(j),omega);
    end
    for j=max(1,dropcount4-3):dropcount4
        % wave field at time tvals(i)=t
        H(((i-1)*length(xvals)+1):(i*length(xvals)),:) = H(((i-1)*length(xvals)+1):(i*length(xvals)),:) ...
            + single_wave(dropxposns4(j),dropyposns4(j),xvals,yvals,t-droptimes4(j),omega);
    end
    for j=max(1,dropcount5-3):dropcount5
        % wave field at time tvals(i)=t
        H(((i-1)*length(xvals)+1):(i*length(xvals)),:) = H(((i-1)*length(xvals)+1):(i*length(xvals)),:) ...
            + single_wave(dropxposns5(j),dropyposns5(j),xvals,yvals,t-droptimes5(j),omega);
    end

    % size(H(i*(1:length(xvals)),:)==0)
    % H(((i-1)*length(xvals)+1):(i*length(xvals)),:) = H(((i-1)*length(xvals)+1):(i*length(xvals)),:) + single_wave(dropxposns(1),dropyposns(1),xvals,yvals,t-droptimes(1),omega);

    % gradient of wave field at this timestep
    [FY,FX] = gradient(H(((i-1)*length(xvals)+1):(i*length(xvals)),:));
    gradH_X(((i-1)*length(xvals)+1):(i*length(xvals)),:) = FX;
    gradH_Y(((i-1)*length(xvals)+1):(i*length(xvals)),:) = FY;
    
    

    % make sure droplet doesn't pass walls
    % if x>=xmax && xprime>0
    %     x=xmax;
    %     xprime=-xprime;
    % elseif x<=xmin && xprime<0
    %     x=xmin;
    %     xprime=-xprime;
    % end
    % if y>=ymax && yprime>0
    %     y=ymax;
    %     yprime=-yprime;
    % elseif y<=ymin && yprime<0
    %     y=ymin;
    %     yprime=-yprime;
    % end


    boundary_radius=xmax-x0;

    % droplet 1
    rad1=sqrt(x1^2+y1^2);
    radprime1=(x1*xprime1+y1*yprime1)/rad1;
    if x1>=0 && y1>=0
        theta1=atan(y1/x1);
    elseif x1<=0
        theta1=atan(y1/x1)+pi;
    elseif x1>=0 && y1<=0
        theta1=atan(y1/x1)+2*pi;
    end
    thetaprime1=(x1*yprime1-y1*xprime1)/(rad1^2);

    %droplet 2
    rad2=sqrt(x2^2+y2^2);
    radprime2=(x2*xprime2+y2*yprime2)/rad2;
    if x2>=0 && y2>=0
        theta2=atan(y2/x2);
    elseif x2<=0
        theta2=atan(y2/x2)+pi;
    elseif x2>=0 && y2<=0
        theta2=atan(y2/x2)+2*pi;
    end
    thetaprime2=(x2*yprime2-y2*xprime2)/(rad2^2);

    %droplet 3
    rad3=sqrt(x3^2+y3^2);
    radprime3=(x3*xprime3+y3*yprime3)/rad3;
    if x3>=0 && y3>=0
        theta3=atan(y3/x3);
    elseif x3<=0
        theta3=atan(y3/x3)+pi;
    elseif x3>=0 && y3<=0
        theta3=atan(y3/x3)+2*pi;
    end
    thetaprime3=(x3*yprime3-y3*xprime3)/(rad3^2);

    rad4=sqrt(x4^2+y4^2);
    radprime4=(x4*xprime4+y4*yprime4)/rad4;
    if x4>=0 && y4>=0
        theta4=atan(y4/x4);
    elseif x4<=0
        theta4=atan(y4/x4)+pi;
    elseif x4>=0 && y4<=0
        theta4=atan(y4/x4)+2*pi;
    end
    thetaprime4=(x4*yprime4-y4*xprime4)/(rad4^2);

    rad5=sqrt(x5^2+y5^2);
    radprime5=(x5*xprime5+y5*yprime5)/rad5;
    if x5>=0 && y5>=0
        theta5=atan(y5/x5);
    elseif x5<=0
        theta5=atan(y5/x5)+pi;
    elseif x5>=0 && y5<=0
        theta5=atan(y5/x5)+2*pi;
    end
    thetaprime5=(x5*yprime5-y5*xprime5)/(rad5^2);

    if rad1>=boundary_radius && radprime1>0 % if going past boundary
        % [rad1,theta1]
        x1=boundary_radius*cos(theta1);
        y1=boundary_radius*sin(theta1);
        radprime1=-radprime1;
        xprime1=rad1*(-sin(theta1))*thetaprime1+radprime1*cos(theta1);
        yprime1=rad1*(cos(theta1))*thetaprime1+radprime1*sin(theta1);
    end
    if rad2>=boundary_radius && radprime2>0 % if going past boundary
        % [rad2,theta2,x2,y2]
        x2=boundary_radius*cos(theta2);
        y2=boundary_radius*sin(theta2);
        radprime2=-radprime2;
        xprime2=rad2*(-sin(theta2))*thetaprime2+radprime2*cos(theta2);
        yprime2=rad2*(cos(theta2))*thetaprime2+radprime2*sin(theta2);
    end
    if rad3>=boundary_radius && radprime3>0 % if going past boundary
        x3=boundary_radius*cos(theta3);
        y3=boundary_radius*sin(theta3);
        radprime3=-radprime3;
        xprime3=rad3*(-sin(theta3))*thetaprime3+radprime3*cos(theta3);
        yprime3=rad3*(cos(theta3))*thetaprime3+radprime3*sin(theta3);
    end
    if rad4>=boundary_radius && radprime4>0 % if going past boundary
        x4=boundary_radius*cos(theta4);
        y4=boundary_radius*sin(theta4);
        radprime4=-radprime4;
        xprime4=rad4*(-sin(theta4))*thetaprime4+radprime4*cos(theta4);
        yprime4=rad4*(cos(theta4))*thetaprime4+radprime4*sin(theta4);
    end
    if rad5>=boundary_radius && radprime5>0 % if going past boundary
        x5=boundary_radius*cos(theta5);
        y5=boundary_radius*sin(theta5);
        radprime5=-radprime5;
        xprime5=rad5*(-sin(theta5))*thetaprime5+radprime5*cos(theta5);
        yprime5=rad5*(cos(theta5))*thetaprime5+radprime5*sin(theta5);
    end


    % solve for horizontal trajectory of drop between this timestep (t) and
    % next timestep (t+dt)
    if i<Nt+1
        if dropcount1>=1
            lastdropposx1 = dropxposns1(dropcount1);
            lastdropposy1 = dropyposns1(dropcount1);
            x_ind1 = find(abs(xvals-lastdropposx1)<=dx/2+0.001,1);
            y_ind1 = find(abs(yvals-lastdropposy1)<=dx/2+0.001,1);
            % x_ind = find(abs(xvals-x)<=dx/2,1);
            % y_ind = find(abs(yvals-y)<=dy/2,1);
            tspan=t:dt/Ntt:t+dt;
            x_init1=[x1; xprime1];
            y_init1=[y1; yprime1];
            [tnew,xnew1] = ode45(@(tt,xx) horizontal_trajectory(tt,xx,gradH_X(i*x_ind1,y_ind1),m,g,D),tspan,x_init1);
            [~,ynew1] = ode45(@(tt,yy) horizontal_trajectory(tt,yy,gradH_Y(i*x_ind1,y_ind1),m,g,D),tspan,y_init1);
            % tnew: Ntt+1 by 1; xnew: Ntt+1 by 2

            full_x1trajectory(1,(i-1)*Ntt+1:i*Ntt) = xnew1(1:Ntt,1);
            full_x1trajectory(2,(i-1)*Ntt+1:i*Ntt) = xnew1(1:Ntt,2);
            % full_x1velocities((i-1)*Ntt+1:i*Ntt) = xnew1(1:Ntt,2);
            full_y1trajectory(1,(i-1)*Ntt+1:i*Ntt) = ynew1(1:Ntt,1);
            full_y1trajectory(2,(i-1)*Ntt+1:i*Ntt) = ynew1(1:Ntt,2);
            % full_y1velocities((i-1)*Ntt+1:i*Ntt) = ynew1(1:Ntt,2);

            nextx1=xnew1(end,1);
            nextxprime1=xnew1(end,2);
            nexty1=ynew1(end,1);
            nextyprime1=ynew1(end,2);
        else
            if xprime1==0
                full_x1trajectory(1,(i-1)*Ntt+1:i*Ntt) = x1;
                full_x1trajectory(2,(i-1)*Ntt+1:i*Ntt) = xprime1;
                % full_x1velocities((i-1)*Ntt+1:i*Ntt) = xprime1;
            else
                full_x1trajectory(1,(i-1)*Ntt+1:i*Ntt) = (x1+xprime1*dt/Ntt):(xprime1*dt/Ntt):(x1+xprime1*dt);
                full_x1trajectory(2,(i-1)*Ntt+1:i*Ntt) = xprime1*ones(Ntt,1);
                % full_x1velocities((i-1)*Ntt+1:i*Ntt) = xprime1*ones(Ntt,1);
            end

            if yprime1==0
                full_y1trajectory(1,(i-1)*Ntt+1:i*Ntt) = y1;
                full_y1trajectory(2,(i-1)*Ntt+1:i*Ntt) = yprime1;
                % full_y1velocities((i-1)*Ntt+1:i*Ntt) = yprime1;
            else
                full_y1trajectory(1,(i-1)*Ntt+1:i*Ntt) = (y1+yprime1*dt/Ntt):(yprime1*dt/Ntt):(y1+yprime1*dt);
                full_y1trajectory(2,(i-1)*Ntt+1:i*Ntt) = yprime1*ones(Ntt,1);
                % full_y1velocities((i-1)*Ntt+1:i*Ntt) = yprime1*ones(Ntt,1);
            end

            nextx1=x1+xprime1*dt;
            nextxprime1=xprime1;
            nexty1=y1+yprime1*dt;
            nextyprime1=yprime1;
        end
        if dropcount2>=1
            lastdropposx2 = dropxposns2(dropcount2);
            lastdropposy2 = dropyposns2(dropcount2);
            x_ind2 = find(abs(xvals-lastdropposx2)<=dx/2+0.001,1);
            y_ind2 = find(abs(yvals-lastdropposy2)<=dx/2+0.001,1);
            % x_ind = find(abs(xvals-x)<=dx/2,1);
            % y_ind = find(abs(yvals-y)<=dy/2,1);
            tspan=t:dt/Ntt:t+dt;
            x_init2=[x2; xprime2];
            y_init2=[y2; yprime2];
            [tnew,xnew2] = ode45(@(tt,xx) horizontal_trajectory(tt,xx,gradH_X(i*x_ind2,y_ind2),m,g,D),tspan,x_init2);
            [~,ynew2] = ode45(@(tt,yy) horizontal_trajectory(tt,yy,gradH_Y(i*x_ind2,y_ind2),m,g,D),tspan,y_init2);
            % tnew: Ntt+1 by 1; xnew: Ntt+1 by 2

            full_x2trajectory(1,(i-1)*Ntt+1:i*Ntt) = xnew2(1:Ntt,1);
            full_x2trajectory(2,(i-1)*Ntt+1:i*Ntt) = xnew2(1:Ntt,2);
            % full_x2velocities((i-1)*Ntt+1:i*Ntt) = xnew2(1:Ntt,2);
            full_y2trajectory(1,(i-1)*Ntt+1:i*Ntt) = ynew2(1:Ntt,1);
            full_y2trajectory(2,(i-1)*Ntt+1:i*Ntt) = ynew2(1:Ntt,2);
            % full_y2velocities((i-1)*Ntt+1:i*Ntt) = ynew2(1:Ntt,2);

            nextx2=xnew2(end,1);
            nextxprime2=xnew2(end,2);
            nexty2=ynew2(end,1);
            nextyprime2=ynew2(end,2);
        else
            if xprime2==0
                full_x2trajectory(1,(i-1)*Ntt+1:i*Ntt) = x2;
                full_x2trajectory(2,(i-1)*Ntt+1:i*Ntt) = xprime2;
                % full_x2velocities((i-1)*Ntt+1:i*Ntt) = xprime2;
            else
                full_x2trajectory(1,(i-1)*Ntt+1:i*Ntt) = (x2+xprime2*dt/Ntt):(xprime2*dt/Ntt):(x2+xprime2*dt);
                full_x2trajectory(2,(i-1)*Ntt+1:i*Ntt) = xprime2*ones(Ntt,1);
                % full_x2velocities((i-1)*Ntt+1:i*Ntt) = xprime2*ones(Ntt,1);
            end

            if yprime2==0
                full_y2trajectory(1,(i-1)*Ntt+1:i*Ntt) = y2;
                full_y2trajectory(2,(i-1)*Ntt+1:i*Ntt) = yprime2;
                % full_y2velocities((i-1)*Ntt+1:i*Ntt) = yprime2;
            else
                full_y2trajectory(1,(i-1)*Ntt+1:i*Ntt) = (y2+yprime2*dt/Ntt):(yprime2*dt/Ntt):(y2+yprime2*dt);
                full_y2trajectory(2,(i-1)*Ntt+1:i*Ntt) = yprime2*ones(Ntt,1);
                % full_y2velocities((i-1)*Ntt+1:i*Ntt) = yprime2*ones(Ntt,1);
            end

            nextx2=x2+xprime2*dt;
            nextxprime2=xprime2;
            nexty2=y2+yprime2*dt;
            nextyprime2=yprime2;
        end
        if dropcount3>=1
            lastdropposx3 = dropxposns3(dropcount3);
            lastdropposy3 = dropyposns3(dropcount3);
            x_ind3 = find(abs(xvals-lastdropposx3)<=dx/2+0.001,1);
            y_ind3 = find(abs(yvals-lastdropposy3)<=dx/2+0.001,1);
            tspan=t:dt/Ntt:t+dt;
            x_init3=[x3; xprime3];
            y_init3=[y3; yprime3];
            [tnew,xnew3] = ode45(@(tt,xx) horizontal_trajectory(tt,xx,gradH_X(i*x_ind3,y_ind3),m,g,D),tspan,x_init3);
            [~,ynew3] = ode45(@(tt,yy) horizontal_trajectory(tt,yy,gradH_Y(i*x_ind3,y_ind3),m,g,D),tspan,y_init3);
            % tnew: Ntt+1 by 1; xnew: Ntt+1 by 2

            full_x3trajectory(1,(i-1)*Ntt+1:i*Ntt) = xnew3(1:Ntt,1);
            full_x3trajectory(2,(i-1)*Ntt+1:i*Ntt) = xnew3(1:Ntt,2);
            full_y3trajectory(1,(i-1)*Ntt+1:i*Ntt) = ynew3(1:Ntt,1);
            full_y3trajectory(2,(i-1)*Ntt+1:i*Ntt) = ynew3(1:Ntt,2);

            nextx3=xnew3(end,1);
            nextxprime3=xnew3(end,2);
            nexty3=ynew3(end,1);
            nextyprime3=ynew3(end,2);
        else
            if xprime3==0
                full_x3trajectory(1,(i-1)*Ntt+1:i*Ntt) = x3;
                full_x3trajectory(2,(i-1)*Ntt+1:i*Ntt) = xprime3;
            else
                full_x3trajectory(1,(i-1)*Ntt+1:i*Ntt) = (x3+xprime3*dt/Ntt):(xprime3*dt/Ntt):(x3+xprime3*dt);
                full_x3trajectory(2,(i-1)*Ntt+1:i*Ntt) = xprime3*ones(Ntt,1);
            end

            if yprime3==0
                full_y3trajectory(1,(i-1)*Ntt+1:i*Ntt) = y3;
                full_y3trajectory(2,(i-1)*Ntt+1:i*Ntt) = yprime3;
            else
                full_y3trajectory(1,(i-1)*Ntt+1:i*Ntt) = (y3+yprime3*dt/Ntt):(yprime3*dt/Ntt):(y3+yprime3*dt);
                full_y3trajectory(2,(i-1)*Ntt+1:i*Ntt) = yprime3*ones(Ntt,1);
            end

            nextx3=x3+xprime3*dt;
            nextxprime3=xprime3;
            nexty3=y3+yprime3*dt;
            nextyprime3=yprime3;
        end
        if dropcount4>=1
            lastdropposx4 = dropxposns4(dropcount4);
            lastdropposy4 = dropyposns4(dropcount4);
            x_ind4 = find(abs(xvals-lastdropposx4)<=dx/2+0.001,1);
            y_ind4 = find(abs(yvals-lastdropposy4)<=dx/2+0.001,1);
            tspan=t:dt/Ntt:t+dt;
            x_init4=[x4; xprime4];
            y_init4=[y4; yprime4];
            [tnew,xnew4] = ode45(@(tt,xx) horizontal_trajectory(tt,xx,gradH_X(i*x_ind4,y_ind4),m,g,D),tspan,x_init4);
            [~,ynew4] = ode45(@(tt,yy) horizontal_trajectory(tt,yy,gradH_Y(i*x_ind4,y_ind4),m,g,D),tspan,y_init4);
            % tnew: Ntt+1 by 1; xnew: Ntt+1 by 2

            full_x4trajectory(1,(i-1)*Ntt+1:i*Ntt) = xnew4(1:Ntt,1);
            full_x4trajectory(2,(i-1)*Ntt+1:i*Ntt) = xnew4(1:Ntt,2);
            full_y4trajectory(1,(i-1)*Ntt+1:i*Ntt) = ynew4(1:Ntt,1);
            full_y4trajectory(2,(i-1)*Ntt+1:i*Ntt) = ynew4(1:Ntt,2);

            nextx4=xnew4(end,1);
            nextxprime4=xnew4(end,2);
            nexty4=ynew4(end,1);
            nextyprime4=ynew4(end,2);
        else
            if xprime4==0
                full_x4trajectory(1,(i-1)*Ntt+1:i*Ntt) = x4;
                full_x4trajectory(2,(i-1)*Ntt+1:i*Ntt) = xprime4;
            else
                full_x4trajectory(1,(i-1)*Ntt+1:i*Ntt) = (x4+xprime4*dt/Ntt):(xprime4*dt/Ntt):(x4+xprime4*dt);
                full_x4trajectory(2,(i-1)*Ntt+1:i*Ntt) = xprime4*ones(Ntt,1);
            end

            if yprime4==0
                full_y4trajectory(1,(i-1)*Ntt+1:i*Ntt) = y4;
                full_y4trajectory(2,(i-1)*Ntt+1:i*Ntt) = yprime4;
            else
                full_y4trajectory(1,(i-1)*Ntt+1:i*Ntt) = (y4+yprime4*dt/Ntt):(yprime4*dt/Ntt):(y4+yprime4*dt);
                full_y4trajectory(2,(i-1)*Ntt+1:i*Ntt) = yprime4*ones(Ntt,1);
            end

            nextx4=x4+xprime4*dt;
            nextxprime4=xprime4;
            nexty4=y4+yprime4*dt;
            nextyprime4=yprime4;
        end
        if dropcount5>=1
            lastdropposx5 = dropxposns5(dropcount5);
            lastdropposy5 = dropyposns5(dropcount5);
            x_ind5 = find(abs(xvals-lastdropposx5)<=dx/2+0.001,1);
            y_ind5 = find(abs(yvals-lastdropposy5)<=dx/2+0.001,1);
            tspan=t:dt/Ntt:t+dt;
            x_init5=[x5; xprime5];
            y_init5=[y5; yprime5];
            [tnew,xnew5] = ode45(@(tt,xx) horizontal_trajectory(tt,xx,gradH_X(i*x_ind5,y_ind5),m,g,D),tspan,x_init5);
            [~,ynew5] = ode45(@(tt,yy) horizontal_trajectory(tt,yy,gradH_Y(i*x_ind5,y_ind5),m,g,D),tspan,y_init5);
            % tnew: Ntt+1 by 1; xnew: Ntt+1 by 2

            full_x5trajectory(1,(i-1)*Ntt+1:i*Ntt) = xnew5(1:Ntt,1);
            full_x5trajectory(2,(i-1)*Ntt+1:i*Ntt) = xnew5(1:Ntt,2);
            full_y5trajectory(1,(i-1)*Ntt+1:i*Ntt) = ynew5(1:Ntt,1);
            full_y5trajectory(2,(i-1)*Ntt+1:i*Ntt) = ynew5(1:Ntt,2);

            nextx5=xnew5(end,1);
            nextxprime5=xnew5(end,2);
            nexty5=ynew5(end,1);
            nextyprime5=ynew5(end,2);
        else
            if xprime5==0
                full_x5trajectory(1,(i-1)*Ntt+1:i*Ntt) = x5;
                full_x5trajectory(2,(i-1)*Ntt+1:i*Ntt) = xprime5;
            else
                full_x5trajectory(1,(i-1)*Ntt+1:i*Ntt) = (x5+xprime5*dt/Ntt):(xprime5*dt/Ntt):(x5+xprime5*dt);
                full_x5trajectory(2,(i-1)*Ntt+1:i*Ntt) = xprime5*ones(Ntt,1);
            end

            if yprime5==0
                full_y5trajectory(1,(i-1)*Ntt+1:i*Ntt) = y5;
                full_y5trajectory(2,(i-1)*Ntt+1:i*Ntt) = yprime5;
            else
                full_y5trajectory(1,(i-1)*Ntt+1:i*Ntt) = (y5+yprime5*dt/Ntt):(yprime5*dt/Ntt):(y5+yprime5*dt);
                full_y5trajectory(2,(i-1)*Ntt+1:i*Ntt) = yprime5*ones(Ntt,1);
            end

            nextx5=x5+xprime5*dt;
            nextxprime5=xprime5;
            nexty5=y5+yprime5*dt;
            nextyprime5=yprime5;
        end
    else
        full_x1trajectory(1,end) = x1;
        full_x1trajectory(2,end) = xprime1;
        full_x2trajectory(1,end) = x2;
        full_x2trajectory(2,end) = xprime2;
        full_x3trajectory(1,end) = x3;
        full_x3trajectory(2,end) = xprime3;
        full_x4trajectory(1,end) = x4;
        full_x4trajectory(2,end) = xprime4;
        full_x5trajectory(1,end) = x5;
        full_x5trajectory(2,end) = xprime5;

        full_y1trajectory(1,end) = y1;
        full_y1trajectory(2,end) = yprime1;
        full_y2trajectory(1,end) = y2;
        full_y2trajectory(2,end) = yprime2;
        full_y3trajectory(1,end) = y3;
        full_y3trajectory(2,end) = yprime3;
        full_y4trajectory(1,end) = y4;
        full_y4trajectory(2,end) = yprime4;
        full_y5trajectory(1,end) = y5;
        full_y5trajectory(2,end) = yprime5;

        nextx1=x1;
        nextxprime1=xprime1;
        nexty1=y1;
        nextyprime1=yprime1;
        nextx2=x2;
        nextxprime2=xprime2;
        nexty2=y2;
        nextyprime2=yprime2;
        nextx3=x3;
        nextxprime3=xprime3;
        nexty3=y3;
        nextyprime3=yprime3;
        nextx4=x4;
        nextxprime4=xprime4;
        nexty4=y4;
        nextyprime4=yprime4;
        nextx5=x5;
        nextxprime5=xprime5;
        nexty5=y5;
        nextyprime5=yprime5;

    end


    % solve for vertical trajectory of drop between this timestep (t) and
    % next timestep (t+dt)
    if i<Nt+1
        x_ind1 = find(abs(xvals-x1)<=dx/2+0.001,1);
        y_ind1 = find(abs(yvals-y1)<=dy/2+0.001,1);
        tspan=t:dt/Ntt:t+dt;
        zb1 = zbar(z1,H(i*x_ind1,y_ind1));
        z_init1 = [z1; zprime1; zb1];
        % calculate time derivative of wave field up to this time
        Hprime1 = gradient(H(x_ind1:length(xvals):i*length(xvals),y_ind1));
        
        [tnew,znew1] = ode45(@(tt,zz) vertical_trajectory(tt,zz,Hprime1(i),m,g,gamma,k,b),tspan,z_init1);

        full_z1trajectory(1,(i-1)*Ntt+1:i*Ntt) = znew1(1:Ntt,1);
        full_z1trajectory(2,(i-1)*Ntt+1:i*Ntt) = znew1(1:Ntt,2);
        % full_z1velocities((i-1)*Ntt+1:i*Ntt) = znew1(1:Ntt,2);

        nextz1=znew1(end,1);
        nextzprime1=znew1(end,2);

        x_ind2 = find(abs(xvals-x2)<=dx/2+0.001,1);
        y_ind2 = find(abs(yvals-y2)<=dy/2+0.001,1);
        tspan=t:dt/Ntt:t+dt;
        zb2 = zbar(z2,H(i*x_ind2,y_ind2));
        z_init2 = [z2; zprime2; zb2];
        % calculate time derivative of wave field up to this time
        Hprime2 = gradient(H(x_ind2:length(xvals):i*length(xvals),y_ind2));
        
        [~,znew2] = ode45(@(tt,zz) vertical_trajectory(tt,zz,Hprime2(i),m,g,gamma,k,b),tspan,z_init2);

        full_z2trajectory(1,(i-1)*Ntt+1:i*Ntt) = znew2(1:Ntt,1);
        full_z2trajectory(2,(i-1)*Ntt+1:i*Ntt) = znew2(1:Ntt,2);
        % full_z2velocities((i-1)*Ntt+1:i*Ntt) = znew2(1:Ntt,2);

        nextz2=znew2(end,1);
        nextzprime2=znew2(end,2);

        x_ind3 = find(abs(xvals-x3)<=dx/2+0.001,1);
        y_ind3 = find(abs(yvals-y3)<=dy/2+0.001,1);
        tspan=t:dt/Ntt:t+dt;
        zb3 = zbar(z3,H(i*x_ind3,y_ind3));
        z_init3 = [z3; zprime3; zb3];
        % calculate time derivative of wave field up to this time
        Hprime3 = gradient(H(x_ind3:length(xvals):i*length(xvals),y_ind3));
        
        [~,znew3] = ode45(@(tt,zz) vertical_trajectory(tt,zz,Hprime3(i),m,g,gamma,k,b),tspan,z_init3);

        full_z3trajectory(1,(i-1)*Ntt+1:i*Ntt) = znew3(1:Ntt,1);
        full_z3trajectory(2,(i-1)*Ntt+1:i*Ntt) = znew3(1:Ntt,2);

        nextz3=znew3(end,1);
        nextzprime3=znew3(end,2);

        x_ind4 = find(abs(xvals-x4)<=dx/2+0.001,1);
        y_ind4 = find(abs(yvals-y4)<=dy/2+0.001,1);
        tspan=t:dt/Ntt:t+dt;
        zb4 = zbar(z4,H(i*x_ind4,y_ind4));
        z_init4 = [z4; zprime4; zb4];
        % calculate time derivative of wave field up to this time
        Hprime4 = gradient(H(x_ind4:length(xvals):i*length(xvals),y_ind4));
        
        [~,znew4] = ode45(@(tt,zz) vertical_trajectory(tt,zz,Hprime4(i),m,g,gamma,k,b),tspan,z_init4);

        full_z4trajectory(1,(i-1)*Ntt+1:i*Ntt) = znew4(1:Ntt,1);
        full_z4trajectory(2,(i-1)*Ntt+1:i*Ntt) = znew4(1:Ntt,2);

        nextz4=znew4(end,1);
        nextzprime4=znew4(end,2);

        x_ind5 = find(abs(xvals-x5)<=dx/2+0.001,1);
        y_ind5 = find(abs(yvals-y5)<=dy/2+0.001,1);
        tspan=t:dt/Ntt:t+dt;
        zb5 = zbar(z5,H(i*x_ind5,y_ind5));
        z_init5 = [z5; zprime5; zb5];
        % calculate time derivative of wave field up to this time
        Hprime5 = gradient(H(x_ind5:length(xvals):i*length(xvals),y_ind5));
        
        [~,znew5] = ode45(@(tt,zz) vertical_trajectory(tt,zz,Hprime5(i),m,g,gamma,k,b),tspan,z_init5);

        full_z5trajectory(1,(i-1)*Ntt+1:i*Ntt) = znew5(1:Ntt,1);
        full_z5trajectory(2,(i-1)*Ntt+1:i*Ntt) = znew5(1:Ntt,2);

        nextz5=znew5(end,1);
        nextzprime5=znew5(end,2);
    else
        full_z1trajectory(1,end) = znew1(end,1);
        full_z1trajectory(2,end) = znew1(end,2);

        nextz1=z1;
        nextzprime1=zprime1;

        full_z2trajectory(1,end) = znew2(end,1);
        full_z2trajectory(2,end) = znew2(end,2);

        nextz2=z2;
        nextzprime2=zprime2;

        full_z3trajectory(1,end) = znew3(end,1);
        full_z3trajectory(2,end) = znew3(end,2);

        nextz3=z3;
        nextzprime3=zprime3;

        full_z4trajectory(1,end) = znew4(end,1);
        full_z4trajectory(2,end) = znew4(end,2);

        nextz4=z4;
        nextzprime4=zprime4;

        full_z5trajectory(1,end) = znew5(end,1);
        full_z5trajectory(2,end) = znew5(end,2);

        nextz5=z5;
        nextzprime5=zprime5;
    end

    
    x_ind1 = find(abs(xvals-x1)<=dx/2+0.001,1);
    y_ind1 = find(abs(yvals-y1)<=dy/2+0.001,1);
    x_ind2 = find(abs(xvals-x2)<=dx/2+0.001,1);
    y_ind2 = find(abs(yvals-y2)<=dy/2+0.001,1);
    x_ind3 = find(abs(xvals-x3)<=dx/2+0.001,1);
    y_ind3 = find(abs(yvals-y3)<=dy/2+0.001,1);
    x_ind4 = find(abs(xvals-x4)<=dx/2+0.001,1);
    y_ind4 = find(abs(yvals-y4)<=dy/2+0.001,1);
    x_ind5 = find(abs(xvals-x5)<=dx/2+0.001,1);
    y_ind5 = find(abs(yvals-y5)<=dy/2+0.001,1);

    % check if drop occurs
    % [z-H(i*x_ind,y_ind),nextz-H(i*x_ind,y_ind)]
    % [z,nextz]
    if z1>=min(H(i*x_ind1,y_ind1),0) && nextz1<=max(H(i*x_ind1,y_ind1),0) && zprime1<=0
        dropxposns1(dropcount1+1)=x1;
        dropyposns1(dropcount1+1)=y1;
        dropzposns1(dropcount1+1)=z1;
        droptimes1(dropcount1+1)=t;
        dropcount1=dropcount1+1;
    end
    if z2>=min(H(i*x_ind2,y_ind2),0) && nextz2<=max(H(i*x_ind2,y_ind2),0) && zprime2<=0
        dropxposns2(dropcount2+1)=x2;
        dropyposns2(dropcount2+1)=y2;
        dropzposns2(dropcount2+1)=z2;
        droptimes2(dropcount2+1)=t;
        dropcount2=dropcount2+1;
    end
    if z3>=min(H(i*x_ind3,y_ind3),0) && nextz3<=max(H(i*x_ind3,y_ind3),0) && zprime3<=0
        dropxposns3(dropcount3+1)=x3;
        dropyposns3(dropcount3+1)=y3;
        dropzposns3(dropcount3+1)=z3;
        droptimes3(dropcount3+1)=t;
        dropcount3=dropcount3+1;
    end
    if z4>=min(H(i*x_ind4,y_ind4),0) && nextz4<=max(H(i*x_ind4,y_ind4),0) && zprime4<=0
        dropxposns4(dropcount4+1)=x4;
        dropyposns4(dropcount4+1)=y4;
        dropzposns4(dropcount4+1)=z4;
        droptimes4(dropcount4+1)=t;
        dropcount4=dropcount4+1;
    end
    if z5>=min(H(i*x_ind5,y_ind5),0) && nextz5<=max(H(i*x_ind5,y_ind5),0) && zprime5<=0
        dropxposns5(dropcount5+1)=x5;
        dropyposns5(dropcount5+1)=y5;
        dropzposns5(dropcount5+1)=z5;
        droptimes5(dropcount5+1)=t;
        dropcount5=dropcount5+1;
    end


    %set new values for next timestep (t+dt)
    x1=nextx1;
    xprime1=nextxprime1;
    y1=nexty1;
    yprime1=nextyprime1;
    z1=nextz1;
    zprime1=nextzprime1;
    x2=nextx2;
    xprime2=nextxprime2;
    y2=nexty2;
    yprime2=nextyprime2;
    z2=nextz2;
    zprime2=nextzprime2;
    x3=nextx3;
    xprime3=nextxprime3;
    y3=nexty3;
    yprime3=nextyprime3;
    z3=nextz3;
    zprime3=nextzprime3;
    x4=nextx4;
    xprime4=nextxprime4;
    y4=nexty4;
    yprime4=nextyprime4;
    z4=nextz4;
    zprime4=nextzprime4;
    x5=nextx5;
    xprime5=nextxprime5;
    y5=nexty5;
    yprime5=nextyprime5;
    z5=nextz5;
    zprime5=nextzprime5;

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
        surf(xplot,yplot,H(((i-1)*length(xvals)+1):((i)*length(xvals)),:),'FaceAlpha',0.5)
        % alpha 0.5
        % view(2)
        hold on
        scatter3([x1,x2,x3,x4,x5],[y1,y2,y3,y4,y5],[z1,z2,z3,z4,z5],'red','filled')
        % scatter3([x1,x2],[y1,y2],[z1,z2],'red','filled')
        shading interp
        colormap('winter')
        cb=colorbar;
        cb.Label.String = 'wave height';
        clim([0,1])
        xlim([xmin,xmax])
        ylim([ymin,ymax])
        zlim([-2,2])
        [X,Y,Z]=cylinder(boundary_radius,50);
        surf(X,Y,Z*4-2,'FaceColor',[0,0,0],'EdgeColor','none','FaceAlpha',0.2)
        % alpha 0.1
        hold off
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Wave Field')
        set(gcf,'color','white')
        if makevid==1
            frame = getframe(wave);
            writeVideo(vidObj,frame);
        end
    end
end
toc

if savematfile==1
    save(strcat(savelocation,int2str(setnum),'_data.mat'),...
        'full_x1trajectory','full_x2trajectory','full_x3trajectory','full_x4trajectory','full_x5trajectory',...
        'full_y1trajectory','full_y2trajectory','full_y3trajectory','full_y4trajectory','full_y5trajectory',...
        'full_z1trajectory','full_z2trajectory','full_z3trajectory','full_z4trajectory','full_z5trajectory')
end

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
% if makeplot==1
%     xfig=figure(2);
%     clf
%     hold on
%     plot(ttvals,full_x1trajectory,'black','LineWidth',2)
%     plot(ttvals,full_x1velocities,'red','LineWidth',2)
%     plot(ttvals,full_x2trajectory,'black','LineWidth',2,'LineStyle','--')
%     plot(ttvals,full_x2velocities,'red','LineWidth',2,'LineStyle','--')
%     hold off
%     xlabel('time')
%     title('X Trajectory')
%     legend('droplet 1 position','droplet 1 velocity','droplet 2 position','droplet 2 velocity')
%     ylabel('x position of droplet')
%     set(gcf,'color','white')
%     if savefigures==1
%         savefig(xfig,strcat(savelocation,'_X',int2str(setnum),'.fig'))
%     end
% end
% plot full y trajectory
% if makeplot==1
%     yfig=figure(3);
%     clf
%     hold on
%     plot(ttvals,full_y1trajectory,'black','LineWidth',2)
%     plot(ttvals,full_y1velocities,'red','LineWidth',2)
%     plot(ttvals,full_y2trajectory,'black','LineWidth',2,'LineStyle','--')
%     plot(ttvals,full_y2velocities,'red','LineWidth',2,'LineStyle','--')
%     hold off
%     xlabel('time')
%     title('Y Trajectory')
%     legend('droplet 1 position','droplet 1 velocity','droplet 2 position','droplet 2 velocity')
%     ylabel('y position of droplet')
%     set(gcf,'color','white')
%     if savefigures==1
%         savefig(yfig,strcat(savelocation,'_Y',int2str(setnum),'.fig'))
%     end
% end
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
    % surface([full_x3trajectory(1,:);full_x3trajectory(1,:)],[full_y3trajectory(1,:);full_y3trajectory(1,:)],...
    %     [zeros(size(full_x3trajectory(1,:)));zeros(size(full_x3trajectory(1,:)))],[speedvals3;speedvals3],...
    %     'facecol','no',...
    %     'edgecol','interp',...
    %     'linew',2);
    % surface([full_x4trajectory(1,:);full_x4trajectory(1,:)],[full_y4trajectory(1,:);full_y4trajectory(1,:)],...
    %     [zeros(size(full_x4trajectory(1,:)));zeros(size(full_x4trajectory(1,:)))],[speedvals4;speedvals4],...
    %     'facecol','no',...
    %     'edgecol','interp',...
    %     'linew',2);
    % surface([full_x5trajectory(1,:);full_x5trajectory(1,:)],[full_y5trajectory(1,:);full_y5trajectory(1,:)],...
    %     [zeros(size(full_x5trajectory(1,:)));zeros(size(full_x5trajectory(1,:)))],[speedvals5;speedvals5],...
    %     'facecol','no',...
    %     'edgecol','interp',...
    %     'linew',2);
    scolors1=speedvals1'*red+(1-speedvals1')*blue;
    scolors2=speedvals2'*red+(1-speedvals2')*blue;
    scolors3=speedvals3'*red+(1-speedvals3')*blue;
    scolors4=speedvals4'*red+(1-speedvals4')*blue;
    scolors5=speedvals5'*red+(1-speedvals5')*blue;
    alphas1=speedvals1;
    alphas1(alphas1<0.5)=alphas1(alphas1<0.5)/5;
    alphas2=speedvals2;
    alphas2(alphas2<0.5)=alphas2(alphas2<0.5)/5;
    alphas3=speedvals3;
    alphas3(alphas3<0.5)=alphas3(alphas3<0.5)/5;
    alphas4=speedvals4;
    alphas4(alphas4<0.5)=alphas4(alphas4<0.5)/5;
    alphas5=speedvals5;
    alphas5(alphas5<0.5)=alphas5(alphas5<0.5)/5;
    scatter(full_x1trajectory(1,:),full_y1trajectory(1,:),[],scolors1,'filled','AlphaData',alphas1,'MarkerFaceAlpha','flat');
    scatter(full_x2trajectory(1,:),full_y2trajectory(1,:),[],scolors2,'filled','AlphaData',alphas2,'MarkerFaceAlpha','flat');
    scatter(full_x3trajectory(1,:),full_y3trajectory(1,:),[],scolors3,'filled','AlphaData',alphas3,'MarkerFaceAlpha','flat')
    scatter(full_x4trajectory(1,:),full_y4trajectory(1,:),[],scolors4,'filled','AlphaData',alphas4,'MarkerFaceAlpha','flat');
    scatter(full_x5trajectory(1,:),full_y5trajectory(1,:),[],scolors5,'filled','AlphaData',alphas5,'MarkerFaceAlpha','flat');
    % s.MarkerFaceAlpha='flat';
    % colormap('cool')
    colormap(flip(redtoblue))
    cb=colorbar;
    % cb.Label.String = 'time';
    % clim([ttvals(1),ttvals(end)])
    cb.Label.String = 'speed';
    clim([0,1])
    % xlim([xmin, xmax])
    % ylim([ymin, ymax])
    scatter3([full_x1trajectory(1,1),full_x2trajectory(1,1),full_x3trajectory(1,1),full_x4trajectory(1,1),full_x5trajectory(1,1)],...
        [full_y1trajectory(1,1),full_y2trajectory(1,1),full_y3trajectory(1,1),full_y4trajectory(1,1),full_y5trajectory(1,1)],...
        zeros(1,5),100,'black','filled')
    scatter3([full_x1trajectory(1,end),full_x2trajectory(1,end),full_x3trajectory(1,end),full_x4trajectory(1,end),full_x5trajectory(1,end)],...
        [full_y1trajectory(1,end),full_y2trajectory(1,end),full_y3trajectory(1,end),full_y4trajectory(1,end),full_y5trajectory(1,end)],...
        zeros(1,5),100,'black','filled','square')
    % plot(full_xtrajectory,full_ytrajectory,'MarkerFaceColor',(ttvals/max(ttvals))'*[1,0,0]+(1-ttvals/max(ttvals))'*[0,0,1],'LineWidth',2)
    % plot(ttvals,full_xvelocities,'red','LineWidth',2)
    th=0:pi/50:2*pi;
    plot(boundary_radius*cos(th),boundary_radius*sin(th),'black')
    hold off
    xlabel('x')
    title('Horizontal Trajectory')
    legend('','','','','','start','end')
    % legend('','','start','end')
    ylabel('y')
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    set(gcf,'color','white')
    if savefigures==1
        savefig(horizontal,strcat(savelocation,'_Horizontal',int2str(setnum),'.fig'))
    end
end

% plot full vertical trajectory
% if makeplot==1
%     vertical=figure(5);
%     clf
%     hold on
%     plot(ttvals,0*ttvals,'black','Linewidth',1)
%     plot(ttvals,full_z1trajectory,'black','LineWidth',2)
%     plot(ttvals,full_z1velocities,'red','LineWidth',2)
%     plot(ttvals,full_z2trajectory,'black','LineWidth',2,'LineStyle','--')
%     plot(ttvals,full_z2velocities,'red','LineWidth',2,'LineStyle','--')
%     hold off
%     xlabel('time')
%     title('Vertical Trajectory')
%     legend('x=0','droplet 1 position','droplet 1 velocity','droplet 2 position', 'droplet 2 velocity')
%     ylabel('vertical position of droplet')
%     set(gcf,'color','white')
%     if savefigures==1
%         savefig(vertical,strcat(savelocation,'_Vertical',int2str(setnum),'.fig'))
%     end
% end

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
