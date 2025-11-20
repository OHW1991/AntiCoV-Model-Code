function distances = calculateDistances(x, y, z, obstacles)
    % Input:
    % x - vector of size (2000, 1) representing x coordinates of particles
    % y - vector of size (2000, 1) representing y coordinates of particles
    % z - vector of size (2000, 1) representing z coordinates of particles
    % obstacles - matrix of size (K, 3) representing x, y, z coordinates of K obstacles
    
    % Number of particles
    numParticles = length(x);
    
    % Number of obstacles
    numObstacles = size(obstacles, 1);
    
    % Reshape particle coordinates into column vectors
    x = x(:);
    y = y(:);
    z = z(:);
    
    % Repeat the particle coordinates for each obstacle
    X_particles = repmat(x, 1, numObstacles);
    Y_particles = repmat(y, 1, numObstacles);
    Z_particles = repmat(z, 1, numObstacles);
    
    % Repeat the obstacle coordinates for each particle
    X_obstacles = repmat(obstacles(:, 1)', numParticles, 1);
    Y_obstacles = repmat(obstacles(:, 2)', numParticles, 1);
    Z_obstacles = repmat(obstacles(:, 3)', numParticles, 1);
    
    % Calculate the squared distances (to avoid unnecessary square root computations)
    distances_squared = (X_particles - X_obstacles).^2 + ...
                        (Y_particles - Y_obstacles).^2 + ...
                        (Z_particles - Z_obstacles).^2;
    
    % Calculate the distances by taking the square root of the squared distances
    distances = sqrt(distances_squared);
end
