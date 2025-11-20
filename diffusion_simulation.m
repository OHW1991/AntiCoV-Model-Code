function diffusion_simulation(N, D, dt, totalTime, wallPos)
    % N         - Number of particles
    % D         - Diffusion coefficient
    % dt        - Time step for the simulation
    % totalTime - Total time to run the simulation
    % wallPos   - Position of the wall (should be 0 or negative for wall at x = 0)
    
    % Initialize particle positions randomly at a certain distance from the wall
    x = abs(normrnd(10, 1, N, 1)); % Initial x position normally distributed around 5 units from the wall
    y = normrnd(0, 1, N, 1); % Initial y position normally distributed
    z = normrnd(0, 1, N, 1); % Initial z position normally distributed
    
    % Calculate the number of steps
    numSteps = totalTime / dt;
    
    % Initialize a counter for particles hitting the wall
    hits = 0;
    
    % Simulation loop
    for t = 1:numSteps
        % Displace particles by a normally distributed step with zero mean and std deviation sqrt(2*D*dt)
        dx = sqrt(2 * D * dt) * randn(N, 1);
        dy = sqrt(2 * D * dt) * randn(N, 1);
        dz = sqrt(2 * D * dt) * randn(N, 1);
        
        x = x + dx;
        y = y + dy;
        z = z + dz;
        
        % Check for particles hitting the wall
        hitIndices = x < wallPos;
        hits = hits + sum(hitIndices);
        
        % Reflect particles that hit the wall
        x(hitIndices) = -x(hitIndices);
    end
    
    % Calculate the rate of hits per unit time
    hitRate = hits / totalTime;
    
    % Display the result
    fprintf('Number of particles hitting the wall per unit time: %.2f\n', hitRate);
end