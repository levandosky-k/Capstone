clear
close all

trajectory_plot=0;
speed_plot=1;
show_startend=1;

path='full_tracks.csv'; %13 drops

Mtable=readtable(path);
M=readmatrix(path);

if strcmp(path,path2)
    times=M(:,2);
end

num_timesteps=7000; % max length(M)
num_droplets=2;
stepsize=1; %1 = don't skip any timesteps

red=[1,0,0];
blue=[0,0,1];
colorLength=length(1:stepsize:num_timesteps);
redtoblue=[linspace(red(1),blue(1),colorLength)',linspace(red(2),blue(2),colorLength)',linspace(red(3),blue(3),colorLength)'];


% trajectory plot
% 
if trajectory_plot==1
figure(1)
clf
hold on;
x=M(1:num_timesteps,:);
for i=1:2:2*num_droplets % each droplet
    if max(x(1:stepsize:num_timesteps,i))<400 && min(x(1:stepsize:num_timesteps,i))>50 %skip weird ones
    scatter(x(1:stepsize:end,i),x(1:stepsize:end,i+1),[],'red','filled','MarkerFaceAlpha',0.1);
    if show_startend==1
    scatter(x(1,i),x(1,i+1),100,'black', 'filled')
    scatter(x(num_timesteps,i),x(num_timesteps,i+1),100,'black', 'filled', 'square')
    end
    end
end
% xlim([0 900]);
% ylim([0 900]);
hold off;
end

if speed_plot==1
% speed plot
s0=zeros(num_timesteps-1,num_droplets);
s=zeros(num_timesteps,num_droplets);
x=M(1:num_timesteps,:); % positions of all droplets from time 1 to num_timesteps
x0=x(1:end-1,:); % positions from time 1 to num_timesteps-1
x1=x(2:end,:); % positions from time 2 to num_timesteps
distsqr=(x1-x0).^2; % squared distance between x coords and between y coords
for i=1:2:2*num_droplets % each droplet
    s0(:,(i+1)/2)=sqrt(distsqr(:,i)+distsqr(:,i+1)); % speed between two points
    % scatter(x(:,i),x(:,i+1),'filled');
end
s(1,:)=s0(1,:);
for j=2:num_timesteps-1
    s(j,:)=(s0(j-1,:)+s0(j,:))/2;
end
s(num_timesteps,:)=s0(num_timesteps-1,:);
snorm=s./(max(s)); % number of timesteps by number of droplets
figure(2)
clf
hold on;
for i=1:2:2*num_droplets % each droplet
    if max(x(1:stepsize:num_timesteps,i))<400 && min(x(1:stepsize:num_timesteps,i))>50
    scolors=snorm(1:stepsize:num_timesteps,(i+1)/2)*red+(1-snorm(1:stepsize:num_timesteps,(i+1)/2))*blue;
    alphas=snorm;
    alphas(alphas<0.5)=alphas(alphas<0.5)/5;
    s=scatter(x(1:stepsize:num_timesteps,i),x(1:stepsize:num_timesteps,i+1),[],scolors,'filled','AlphaData',alphas(:,(i+1)/2));
    % s=scatter(x(1:stepsize:num_timesteps,i),x(1:stepsize:num_timesteps,i+1),[],scolors,'filled');
    s.MarkerFaceAlpha='flat';
    xlabel('x position')
    ylabel('y position')
    colormap(flip(redtoblue))
    c = colorbar;
    c.Label.String = 'Normalized Speed';
    title('Trajectory and Speed of Droplet')
    if show_startend==1
    scatter(x(1,i),x(1,i+1),100,'black', 'filled')
    scatter(x(num_timesteps,i),x(num_timesteps,i+1),100,'black', 'filled', 'square')
    legend('','start', 'end')
    end
    end
end
hold off;
end
