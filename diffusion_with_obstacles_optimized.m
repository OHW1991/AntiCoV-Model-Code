function [hitRate] = diffusion_with_obstacles_optimized(N, D, dt, totalTime, wallPos, K, obstacleSize, attachmentTime)
    % N             - Number of particles
    % D             - Diffusion coefficient
    % dt            - Time step for the simulation
    % totalTime     - Total time to run the simulation
    % wallPos       - Position of the wall (should be 0 or negative for wall at x = 0)
    % K             - Number of obstacles
    % obstacleSize  - Size of the obstacles (radius)
    % attachmentTime - Time steps a particle remains attached to an obstacle
    
    % Initialize particle positions randomly at a certain distance from the wall
    x = abs(normrnd(10, 1, N, 1)); % Initial x position normally distributed around 5 units from the wall
    y = normrnd(0, 1, N, 1); % Initial y position normally distributed
    z = normrnd(0, 1, N, 1); % Initial z position normally distributed
    
    % Initialize obstacle positions randomly in the space
    obstaclePositions = normrnd(5, 2.5, K, 3); % Obstacle positions normally distributed in space

    % Initialize particle attachment status
    attached = false(N, 1);
    attachmentCounter = zeros(N, 1);
    
    % Calculate the number of steps
    numSteps = totalTime / dt;
    
    % Initialize a counter for particles hitting the wall
    hits = 0;
    
    % Pre-compute sqrt(2 * D * dt)
    sqrt2Ddt = sqrt(2 * D * dt);
    
    % Simulation loop
    for t = 1:numSteps
        % Displace particles by a normally distributed step with zero mean and std deviation sqrt(2*D*dt)
        dx = sqrt2Ddt * randn(N, 1);
        dy = sqrt2Ddt * randn(N, 1);
        dz = sqrt2Ddt * randn(N, 1);
        
        % Update positions for non-attached particles
        nonAttached = ~attached;
        x(nonAttached) = x(nonAttached) + dx(nonAttached);
        y(nonAttached) = y(nonAttached) + dy(nonAttached);
        z(nonAttached) = z(nonAttached) + dz(nonAttached);

        % Check for particles hitting the wall
        hitIndices = x < wallPos & nonAttached;
        hits = hits + nnz(hitIndices);
        
        % Reflect particles that hit the wall
        x(hitIndices) = -x(hitIndices);
        
        % Check for particles hitting obstacles
        for j = 1:K
            dist = sqrt((x - obstaclePositions(j, 1)).^2 + ...
                        (y - obstaclePositions(j, 2)).^2 + ...
                        (z - obstaclePositions(j, 3)).^2);
            hitObstacles = dist < obstacleSize & nonAttached;
            attached(hitObstacles) = true;
            attachmentCounter(hitObstacles) = attachmentTime;
        end
        
        % Update attachment status
        attached(attachmentCounter > 0) = true;
        attached(attachmentCounter <= 0) = false;
        attachmentCounter(attached) = attachmentCounter(attached) - 1;
        
        % Release particles after attachment time
        releaseIndices = attachmentCounter == 0 & attached;
        if any(releaseIndices)
            % Release particles in a random direction
            theta = 2 * pi * rand(sum(releaseIndices), 1);
            phi = acos(2 * rand(sum(releaseIndices), 1) - 1);
            x(releaseIndices) = x(releaseIndices) + obstacleSize * sin(phi) .* cos(theta);
            y(releaseIndices) = y(releaseIndices) + obstacleSize * sin(phi) .* sin(theta);
            z(releaseIndices) = z(releaseIndices) + obstacleSize * cos(phi);
        end
    end
    
    % Calculate the rate of hits per unit time
    hitRate = hits / totalTime;
    
%     % Display the result
%     fprintf('Number of particles hitting the wall per unit time: %.2f\n', hitRate);
end
