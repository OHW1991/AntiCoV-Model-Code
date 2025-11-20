function [hitRate,particleRtrajectory] = diffusion_with_obstacles_finite_space_optimized(N, D, dt, totalTime, wallPos, K, obstacleSize, attachmentTime, spaceSize)
    % N             - Number of particles
    % D             - Diffusion coefficient
    % dt            - Time step for the simulation
    % totalTime     - Total time to run the simulation
    % wallPos       - Position of the wall (should be 0 for wall at x = 0)
    % K             - Number of obstacles
    % obstacleSize  - Size of the obstacles (radius)
    % attachmentTime - Time steps a particle remains attached to an obstacle
    % spaceSize     - Size of the finite space [x, y, z]
    
    % Initialize particle positions randomly within the finite space
    x = abs(normrnd(spaceSize(1)/2, spaceSize(1)/5, N, 1)); % Initial x position normally distributed
    y = normrnd(spaceSize(2)/2, spaceSize(2)/5, N, 1); % Initial y position normally distributed
    z = normrnd(spaceSize(3)/2, spaceSize(3)/5, N, 1); % Initial z position normally distributed
    
    % Ensure particles start within bounds
    x = min(max(x, 0), spaceSize(1));
    y = min(max(y, 0), spaceSize(2));
    z = min(max(z, 0), spaceSize(3));
    
    
    % Initialize obstacle positions randomly within the space
    obstaclePositions = zeros(K,3); 
    % Obstacle positions uniformly distributed across from the flux wall
    % (x=0)
%     obstaclePositions(:,1) = 0.25 * spaceSize(1) + 0.5 * spaceSize(1) .* rand(K, 1); 
    obstaclePositions(:,1) = spaceSize(1) .* rand(K, 1); 
    obstaclePositions(:,2:3) = spaceSize(2:3) .* rand(K, 2); 
    
    % Ensure obstacle positions are within bounds
    obstaclePositions(:,1) = min(max(obstaclePositions(:, 1)-obstacleSize, 0), spaceSize(1)-obstacleSize);
    obstaclePositions(:,2) = min(max(obstaclePositions(:, 2)-obstacleSize, 0), spaceSize(2)-obstacleSize);
    obstaclePositions(:,3) = min(max(obstaclePositions(:, 3)-obstacleSize, 0), spaceSize(3)-obstacleSize);

    % Initialize particle attachment status
    attached = false(N, 1);
    attachmentCounter = zeros(N, 1);
    
    % Calculate the number of steps
    numSteps = totalTime / dt;
    
    % Initialize matrix for storing particle trajectories
    particleRtrajectory = zeros(numSteps+1,N,3);
    particleRtrajectory(1,:,:) = [x,y,z];
    
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

        % Reflect particles that hit the walls of the finite space
        x(x > spaceSize(1)) = 2 * spaceSize(1) - x(x > spaceSize(1));
        y(y < 0) = -y(y < 0);
        y(y > spaceSize(2)) = 2 * spaceSize(2) - y(y > spaceSize(2));
        z(z < 0) = -z(z < 0);
        z(z > spaceSize(3)) = 2 * spaceSize(3) - z(z > spaceSize(3));

        % Check for particles hitting the wall at wallPos
        hitIndices = x < wallPos & nonAttached;
        hits = hits + nnz(hitIndices);
        
        % Teleport particles that hit the measurement wall (to keep
        % the flux constant)
        x(hitIndices) = abs(normrnd(spaceSize(1)/2, spaceSize(1)/5, sum(hitIndices), 1));
        
        % Check for particles hitting obstacles
        dist = calculateDistances(x,y,z,obstaclePositions);
        hitObstacles = min(dist,[],2) < obstacleSize & nonAttached;
        attached(hitObstacles) = true;
        attachmentCounter(hitObstacles) = attachmentTime;
%         for j = 1:K
%             dist = sqrt((x - obstaclePositions(j, 1)).^2 + ...
%                         (y - obstaclePositions(j, 2)).^2 + ...
%                         (z - obstaclePositions(j, 3)).^2);
%             hitObstacles = dist < obstacleSize & nonAttached;
%             attached(hitObstacles) = true;
%             attachmentCounter(hitObstacles) = attachmentTime;
%         end
        
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
            
            % Ensure released particles are within bounds
            x = min(max(x, 0), spaceSize(1));
            y = min(max(y, 0), spaceSize(2));
            z = min(max(z, 0), spaceSize(3));
            
            
        end
%         scatter3(x,y,z); hold on;
%         scatter3(obstaclePositions(:,1),obstaclePositions(:,2),obstaclePositions(:,3),10,'filled');
%         ph=patch(1e4.*[1 1 1 1],1e4.*[1 0 0 1],1e4.*[1 1 0 0],0.5*1e4.*[1 1 1 1],'facecolor','g','FaceAlpha',0.25,'edgecolor','none');
%         ph=patch(1e4.*[1 0 0 1],1e4.*[1 1 1 1],1e4.*[1 1 0 0],0.5*1e4.*[1 1 1 1],'facecolor','g','FaceAlpha',0.25,'edgecolor','none');
%         ph=patch(1e4.*[1 1 0 0],1e4.*[1 0 0 1],1e4.*[0 0 0 0],0.5*1e4.*[1 1 1 1],'facecolor','b','FaceAlpha',0.25,'edgecolor','none');
%         ph=patch(1e4.*[1 1 0 0],1e4.*[1 0 0 1],1e4.*[1 1 1 1],0.5*1e4.*[1 1 1 1],'facecolor','b','FaceAlpha',0.25,'edgecolor','none');
%         view(3);
%         hold off;
%         xlim([0,1e4]); ylim([0,1e4]); zlim([0,1e4]);
%         xlabel('X');
%         ylabel('Y');
%         zlabel('Z');
%         fprintf('something\n');
    particleRtrajectory(t+1,:,:) = [x,y,z];
    end
    
    % Calculate the rate of hits per unit time
    hitRate = hits / totalTime;
    
end
