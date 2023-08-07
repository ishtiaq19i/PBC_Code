%% Particle Track Code

% x,y,z,u1,v1,w1,x0,y0,z0,n,r,delta_t,method_interp are
% needed as input
% x - grid points in x
% y - grid points in y
% z - grid points in z
% u1 - velocity in x
% v1 - velocity in y
% w1 - velocity in z
% x0, y0, z0 - initial points
% n - number of particles
% r - the number of time steps
% delta_t - interval between consecutive timesteps
% method_interp - interpolation method used to interpolate velocities at
% non-nodal points where velocity data is unavailable

function[x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory] = Particle_Track_Core_RK2_PBC(x,y,z,u1_t1,u1_t2,v1_t1,v1_t2,w1_t1,w1_t2,x0,y0,z0,n,delta_t,x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory)

X = x;
Y = y;
Z = z;

Lx = max(x);
Ly = max(y);
Lz = max(z);

parfor j=1:n

    u_interp1 = interp3(X,Y,Z,u1_t1,x0(j),y0(j),z0(j));
    v_interp1 = interp3(X,Y,Z,v1_t1,x0(j),y0(j),z0(j));
    w_interp1 = interp3(X,Y,Z,w1_t1,x0(j),y0(j),z0(j));

    u_trajectory1 = u_interp1;
	v_trajectory1 = v_interp1;
    w_trajectory1 = w_interp1;
    
	x_track = x0(j) + u_interp1*delta_t;
    y_track = y0(j) + v_interp1*delta_t;
    z_track = z0(j) + w_interp1*delta_t;
    
    % Periodic Boundary Conditions

    if x_track > max(x) || x_track < min(x) || y_track > max(y) || y_track < min(x) || z_track > max(z) || z_track < min(x)
        
        if x_track >= max(x)
            x_track = x_track - (Lx-min(x));
            x0(j) = x0(j) - (Lx-min(x));
        elseif x_track <= min(x)
            x_track = x_track + (Lx-min(x));
            x0(j) = x0(j) + (Lx-min(x));
        end
 
        if y_track >= max(y)
            y_track = y_track - (Ly-min(y));
            y0(j) = y0(j) - (Ly-min(y));
        elseif y_track <= min(y)
            y_track = y_track + (Ly-min(y));
            y0(j) = y0(j) + (Ly-min(y));
        end
 
        if z_track >= max(z)
            z_track = z_track - (Lz-min(z));
            z0(j) = z0(j) - (Lz-min(z));
        elseif z_track <= min(z)
            z_track = z_track + (Lz-min(z));
            z0(j) = z0(j) + (Lz-min(z));
        end
        else
    end
        
    x_trajectory_temp = x_track; % saving trajectory values to prevent zero-ing of data due to using parfor loop
    y_trajectory_temp = y_track;
    z_trajectory_temp = z_track;
    
    u_interp = interp3(X,Y,Z,u1_t2,x_trajectory_temp,y_trajectory_temp,z_trajectory_temp);
    v_interp = interp3(X,Y,Z,v1_t2,x_trajectory_temp,y_trajectory_temp,z_trajectory_temp);
    w_interp = interp3(X,Y,Z,w1_t2,x_trajectory_temp,y_trajectory_temp,z_trajectory_temp);
   
    x_trajectory(j) = x0(j) + (u_trajectory1+u_interp)*delta_t/2;
    y_trajectory(j) = y0(j) + (v_trajectory1+v_interp)*delta_t/2;
    z_trajectory(j) = z0(j) + (w_trajectory1+w_interp)*delta_t/2;
        
    u_trajectory(j) = u_trajectory1;
    v_trajectory(j) = v_trajectory1;
    w_trajectory(j) = w_trajectory1;

end

end




