Lx = 1; % set length in x of box 
        % (this instance is normalized to 1 m)
Ly = 1; % set length in y of box 
        % (this instance is normalized to 1)

Nx = 100; % Nx is the number of intervals that we will divide Lx by
Ny = 100; % Ny is the number of intervals that we will divide Ly by

nx = Nx + 1; % set the number of indices in x (101, from 0 to 100) 
ny = Ny + 1; % set the number of indices in y (101, from 0 to 100)

dx = Lx/Nx; % the physical size of each interval between any
            % x(i)th to x(i+1)th point, etc. (0.01 m)
            
dy = Ly/Ny; % the physical size of each interval between any
            % y(i)th to y(i+1)th point, etc. (0.01 m)

x = (0:Nx)*dx; % the physical values of length along the x-axis
               % (0 m, 0.01 m, 0.02 m, ..., 0.98 m, 0.99 m, 1 m)
y = (0:Ny)*dy; % the physical values of length along the y-axis
               % (0 m, 0.01 m, 0.02 m, ..., 0.98 m, 0.99 m, 1 m)

% We can use the below definitons of i and j 
% to select all points on the grid at once.

i = 2:nx-1; % set a list of indices from the second to the second-last
            % value in the x range (point 1 to 100 along x)
            
j = 2:ny-1; % set a list of indices from the second to the second-last
            % value in the y range (point 1 to 100 along y)
            

V = zeros(nx,ny); % define V as a matrix of zeros of size 101*101

V_old = V; % differentiate the old and new (default) iteration of V 

% Set boundary conditions for V

% V(:,1) = 0; % y min bound
% V(:,ny) = 0; % y max bound

V(1,:) = 1; % x min bound
V(nx,:) = 0; % x max bound

eps = 1.e-6 ; % set eps 

error = 2*eps; % set error criteria
n = 0;

while (error > eps)
    n = n + 1; % set count

    V(i,j) = ( V(i+1,j) ...
             + V(i-1,j) ...
             + V(i,j+1) ...
             + V(i,j-1) )./4 ; % Applying the del^2 operator 
                               % in FD form
                               
    V(:,1) = V(:,2); % for the x min bound boundary, do this to remove BC
    V(:,ny) = V(:,Ny-1); % for y min bound right boundary, do this to remove BC
    
    % the above two lines just fixes the left/right bounds to be the
    % same value as the next left/right values respectively. It effectively
    % removes the zero boundary conditions
    
%     V(:,1) = V(:,2); % for the y min bound, do this to remove BC
%     V(:,ny) = V(:,Ny-1); % for the y max bound, do this to remove BC
                               
    error = max(abs(V(:) - V_old(:)));
     
%      if any(isnan(V(:))) || any(isinf(V(:)))
%          fprintf('iterations diverge \n');
%          return;
%      end
     
     V_old = V ; % refresh
     % fprintf ('%g %e', n, error);
     
     % GRADIENT STUFF
     [Ex,Ey] = gradient(V_old);
     Ex = Ex .*-1;
     Ey = Ey .*-1;
     
     vBox = imboxfilt(V, 3);
     
     % PLOT GAMING
     if mod (n, 10) == 0
        
        subplot(1,3,1)
        surf(V_old)
        
        subplot(1,3,2)
        quiver(Ex, Ey)
        
        subplot(1,3,3)
        imshowpair(V, vBox)
        
        pause(0.05)
     end
     
     % PLOT GAMING

end
fprintf('%g\n', n);
