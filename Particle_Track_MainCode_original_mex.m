tic % measuring the code's elapsed time; starting time

clear

% Please include ReadDNSData3, Particle_Trake_Intermediate,
% Particle_Track_Core_RK2_032923_mex together with this script in the same
% folder. Please note that the directories need to adjust accordingly to
% the location of the data and mex function


%% Importing VOF timestep 1
% to make sure coordinates inputted don't overlap with the bubbles
myfilename = sprintf('Matlab_stag001.dat');
directory = 'G:\Bubbles_DNS\2_05'; % Path can be changed

filename = fullfile(directory, myfilename); % LOCATION OF THE DATA IS CHANGED HERE
[Data] = ReadDNSData3(filename);

vof_1(:,:,:) = Data.vof;

%% Input Values

clear vars x y z 
delta_t = 0.001;

n = 10000; % NUMBER OF PARTICLES
n = single(n);
number_of_timesteps = 542; % NUMBER OF TIMESTEPS

x = double((1:1:330));
y = double((1:1:165));
z = double((1:1:165));
[X,Y,Z] = meshgrid(x,y,z);

x = transpose(x);
y = transpose(y);
z = transpose(z);

mat = single((1:1:165*330*165));
mat_randIdcs = randperm(length(mat),n);
X = reshape(X,[8984250,1]);
Y = reshape(Y,[8984250,1]);
Z = reshape(Z,[8984250,1]);

%generate the random indices

x_randIdcs = randperm(length(X),n);
y_randIdcs = randperm(length(Y),n);
z_randIdcs = randperm(length(Z),n);

x0 = X(x_randIdcs);
y0 = Y(y_randIdcs);
z0 = Z(z_randIdcs);

coordinates_initial = [x0 y0 z0];

parfor i = 1:numel(x0)
    if vof_1(y0(i),x0(i),z0(i)) == 1
        disp(i)
         x_randIdcs = randperm(length(X),1);
         y_randIdcs = randperm(length(Y),1);
         z_randIdcs = randperm(length(Z),1);
         x0(i) = X(x_randIdcs);
         y0(i) = Y(y_randIdcs);
         z0(i) = Z(z_randIdcs);
    else
    end
end

clear vof_1

x0 = x0*(40/165); %from indices to real world coordinates in mm
y0 = y0*(40/165);
z0 = z0*(40/165);

% x_decimal = x_decimal*(40/165); %from indices to real world coordinates in mm
% y_decimal = y_decimal*(40/165);
% z_decimal = z_decimal*(40/165);
%using index values

x0_original = x0;
y0_original = y0;
z0_original = z0;

x = x*(40/165); %from indices to real world coordinates in mm
y = y*(40/165);
z = z*(40/165);

clear x_randIdcs y_randIdcs z_randIdcs

%% Reading Data

q = 1;

x_track = nan(n,1); %n = number of particles; row
y_track = nan(n,1);
z_track = nan(n,1);

u_track = nan(n,1); %n = number of particles; row
v_track = nan(n,1);
w_track = nan(n,1);

%GLOBAL INITIALISATION 

x_trajectory = nan(n,1); %n = number of particles; row
y_trajectory = nan(n,1);
z_trajectory = nan(n,1);

u_trajectory = nan(n,1);
v_trajectory = nan(n,1);
w_trajectory = nan(n,1);

u1 = zeros(165,330,165,2);
v1 = zeros(165,330,165,2);
w1 = zeros(165,330,165,2);
vof = zeros(165,330,165,2);

temp_u1 = zeros(size(u1));
temp_v1 = zeros(size(temp_u1));
temp_w1 = zeros(size(temp_u1));
temp_vof = zeros(size(temp_u1));

% Main Loop
for b = 1:number_of_timesteps

if b == number_of_timesteps
        break
else
r = q + 1;

for a=q:r

    if a <= 9
        myfilename1 = sprintf('Matlab_stag00%d.dat', a);
        
    elseif a > 9 && a <= 99
        myfilename1 = sprintf('Matlab_stag0%d.dat', a);
    else
        myfilename1 = sprintf('Matlab_stag%d.dat', a);
        
    end
    filename1 = fullfile(directory, myfilename1); % LOCATION OF THE DATA IS CHANGED HERE


[Data] = ReadDNSData3(filename1);

if a == q
    temp_u1(:,:,:,1) = Data.u;
    temp_v1(:,:,:,1) = Data.v;
    temp_w1(:,:,:,1) = Data.w;
    temp_vof(:,:,:,1) = Data.vof;
else
    %if a == r
    temp_u1(:,:,:,2) = Data.u;
    temp_v1(:,:,:,2) = Data.v;
    temp_w1(:,:,:,2) = Data.w;
    temp_vof(:,:,:,2) = Data.vof;
end
end

u1(:,:,:,1) = temp_u1(:,:,:,1);
u1(:,:,:,2) = temp_u1(:,:,:,2);
u1 = u1*1000;

v1(:,:,:,1) = temp_v1(:,:,:,1);
v1(:,:,:,2) = temp_v1(:,:,:,2);
v1 = v1*1000;

w1(:,:,:,1) = temp_w1(:,:,:,1);
w1(:,:,:,2) = temp_w1(:,:,:,2);
w1 = w1*1000;

vof(:,:,:,1) = temp_vof(:,:,:,1);
vof(:,:,:,2) = temp_vof(:,:,:,2);

%velocity data is given in m/s
%convert from m/s to mm/s

clear vars temp_u1 temp_v1 temp_w1 temp_vof

%[x_track_int, y_track_int, z_track_int, u_track_int, v_track_int, w_track_int] = Particle_Track_Intermediate(x,y,z,u1,v1,w1,x0,y0,z0,n,delta_t,x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,vof);

[x_track_int, y_track_int, z_track_int, u_track_int, v_track_int, w_track_int] = Particle_Track_Intermediate_PBC(x,y,z,u1,v1,w1,x0,y0,z0,n,delta_t,x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,vof);

q = q + 1;

x0 = x_track_int;
y0 = y_track_int;
z0 = z_track_int;

x_track(:,b) = x_track_int;
y_track(:,b) = y_track_int;
z_track(:,b) = z_track_int;
u_track(:,b) = u_track_int;
v_track(:,b) = v_track_int;
w_track(:,b) = w_track_int;

disp(['Timestep(s) Completed: ' num2str(b)]);

% clearvars -except r q b x y z u1 v1 w1 vof x0_original y0_original z0_original x_track y_track z_track x0 y0 z0 n delta_t
clear vars u1 v1 w1 vof
end % repitition ends
end

x_trajectory = [x0_original x_track];
y_trajectory = [y0_original y_track];
z_trajectory = [z0_original z_track];

u_trajectory = u_track;
v_trajectory = v_track;
w_trajectory = w_track;

clear vars x_track y_track z_track u_track v_track w_track x_track_int y_track_int z_track_int u_track_int v_track_int w_track_int
%%
% for i = 1:number_of_timesteps
%     centroid_vof{i} = regionprops3(vof(:,:,:,i),'Centroid');
% end
%% Plotting Particle Trajectories in 3D
% 
% 
x_track = x_trajectory;
y_track = y_trajectory;
z_track = z_trajectory;

%trajectory = [x_track; y_track; z_track];

F1 = figure;
for i=1:n
    plot3(x_track(i,:),y_track(i,:),z_track(i,:))
    %hold on
    %plot3(x_track_GP(i,:),y_track_GP(i,:),z_track_GP(i,:))
    %hold on
    %plot3(x_track_RK2(i,:),y_track_RK2(i,:),z_track_RK2(i,:))
    axis square
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    %L=legend('particle track','GetPosition','particle track with RK2','location','northeast')

end

%     hold on
% xlim([0 max(x_trajectory)])
% ylim([0 max(y_trajectory)])
% zlim([0 max(z_trajectory)])

hold on 


% p1 = patch(isosurface(X,Y,Z,vof(:,:,:,1),0.5));  
% 
% %isonormals(x,y,z,v,p)
% 
% p1.FaceColor = 'red';
% 
% p1.EdgeColor = 'none';
% 
% daspect([1 1 1])
% 
% view(3); 
% 
% axis tight, grid on 
% 
% camlight

% hold on

% quiver3(X,Y,Z,u(:,:,:,1),v(:,:,:,1),w(:,:,:,1))

% grid on

    xlabel('x (mm)')

    ylabel('y (mm)')

    zlabel('z (mm)')
    


%     

%  contourfm(u)

% clim([min(min(u(:,:,:,1))) max(max(u(:,:,:,1)))])

% cb = contourcbar("eastoutside");
 
% cb.XLabel.String = "Velocity (mm/s)";

% cb.Ticks = levels;

 

toc

