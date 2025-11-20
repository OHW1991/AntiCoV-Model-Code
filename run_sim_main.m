% Run simulation without obstacles 20 times to get statistics
flux_0_obs = zeros(10,1);
for i = 1:10
    flux_0_obs(i,1) = diffusion_with_obstacles_finite_space_optimized(2000, 1e6, 0.01, 10.0, 0, 0, 100, 24, [1e4,1e4,1e4]);
end

fprintf('10 particles sim\n');
flux_10_obs = zeros(10,1);
for i = 1:10
    flux_10_obs(i,1) = diffusion_with_obstacles_finite_space_optimized(2000, 1e6, 0.01, 10.0, 0, 10, 100, 24, [1e4,1e4,1e4]);
end

fprintf('100 particles sim\n');
flux_100_obs = zeros(10,1);
for i = 1:10
    flux_100_obs(i,1) = diffusion_with_obstacles_finite_space_optimized(2000, 1e6, 0.01, 10.0, 0, 100, 100, 24, [1e4,1e4,1e4]);
end

fprintf('500 particles sim\n');
flux_500_obs = zeros(10,1);
for i = 1:10
    flux_500_obs(i,1) = diffusion_with_obstacles_finite_space_optimized(2000, 1e6, 0.01, 10.0, 0, 500, 100, 24, [1e4,1e4,1e4]);
end

fprintf('1000 particles sim\n');
flux_1000_obs = zeros(10,1);
for i = 1:10
    flux_1000_obs(i,1) = diffusion_with_obstacles_finite_space_optimized(2000, 1e6, 0.01, 10.0, 0, 1000, 100, 24, [1e4,1e4,1e4]);
end


boxplot([flux_0_obs,flux_10_obs,flux_100_obs,flux_500_obs,flux_1000_obs],'Labels',{'0','10','100','500','1000'});
xlabel('Number of obstacles');
ylabel('Flux');


 
% Generate heatmap
obstacle_numbers = floor([2e3,1e3,1e3/(3^1),1e3/(3^2),1e3/(3^3),1e3/(3^4),1e3/(3^5),...
    1e3/(3^6),1e3/(3^7)]);
sialic_acid_penalty = 3*[0:4:24];

mean_matrix = zeros(numel(obstacle_numbers),numel(sialic_acid_penalty));
std_matrix = zeros(numel(obstacle_numbers),numel(sialic_acid_penalty));
flux_matrix = zeros(10,numel(obstacle_numbers),numel(sialic_acid_penalty));
for num = 1:numel(obstacle_numbers) 
   for penalty = 1:numel(sialic_acid_penalty)
%        flux_vector = zeros(10,1);
       for iteration = 1:10
%            flux_vector(iteration,1) = diffusion_with_obstacles_finite_space_optimized(2000, 1e6, 0.01, 10.0, 0, obstacle_numbers(num), 100, sialic_acid_penalty(penalty), [1e4,1e4,1e4]);
           flux_matrix(iteration,num,penalty) = diffusion_with_obstacles_finite_space_optimized(2000, 1e6, 0.01, 10.0, 0, obstacle_numbers(num), 100, sialic_acid_penalty(penalty), [1e4,1e4,1e4]);
       end 
       fprintf('iteration over\n');
%        mean_matrix(num,penalty) = mean(flux_vector);
%        std_matrix(num,penalty) = std(flux_vector);
   end
end

heatmap(obstacle_numbers,sialic_acid_penalty/3,mean_matrix');
% caxis([10,35]);
colormap turbo;
xlabel('# Obstacles');
ylabel('# Sialic Acids');

mean_matrix2 = mean_matrix;
hHM= heatmap(obstacle_numbers,sialic_acid_penalty/3,mean_matrix2'/max(mean_matrix2(:)));
caxis([0.5,1]);
colormap turbo;
xlabel('# Obstacles');
ylabel('# Sialic Acids');
hHM.YDisplayData=flip(hHM.YDisplayData);  
hHM.XDisplayData=flip(hHM.XDisplayData);  


% plot boxplots by row
for num=1:numel(sialic_acid_penalty)
   figure;
   boxplot([flux_matrix(:,9,num),flux_matrix(:,8,num),flux_matrix(:,7,num),...
       flux_matrix(:,6,num),flux_matrix(:,5,num),flux_matrix(:,4,num),...
       flux_matrix(:,3,num),flux_matrix(:,2,num),flux_matrix(:,1,num)],...
       'Labels',{'0','1','4','12','37','111','333','1000','2000'});
   xlabel('Number of obstacles');
   ylabel('Flux');
   ylim([40,90]);
   title([num2str(sialic_acid_penalty(num)/3) ' sialic acids']);
end




% plot trajectories
[hhh,trajectories_mat] = diffusion_with_obstacles_finite_space_optimized(2000, 1e6, 0.01, 10.0, 0, 2000, 100, 72, [1e4,1e4,1e4]);
time_vector = 0:1000;

for i = 1:40
    traj_num = randi(2000);
    if (~isempty(find(diff(trajectories_mat(:,traj_num,1))>1000,1)))
        cutoff_point = find(diff(trajectories_mat(:,traj_num,1))>1000,1);
        figure;
        subplot(3,1,1);
        plot(time_vector(1:cutoff_point),trajectories_mat(1:cutoff_point,traj_num,1));
        title('Particle',traj_num);
        ylim([0,1e4]); ylabel('X position'); xlabel('Time Steps');
        subplot(3,1,2);
        plot(time_vector(1:cutoff_point),trajectories_mat(1:cutoff_point,traj_num,2));
        ylim([0,1e4]); ylabel('Y position'); xlabel('Time Steps');
        subplot(3,1,3);
        plot(time_vector(1:cutoff_point),trajectories_mat(1:cutoff_point,traj_num,3));
        ylim([0,1e4]);
        xlabel('Time Steps');
        ylabel('Z position');
    end
end


tic
penalty = 50;
num_runs = 10;
new_particles = 20;
bounding_box = [1e4,1e4,1e4];
D = (0.01*1e4)^2/(2*0.01);
flux_0_obs = zeros(num_runs,1);
particle_counter_0 = zeros(1000,num_runs); 
infection_counter_0_obs = zeros(num_runs,1);
death_counter_0_obs = zeros(num_runs,1);
for i = 1:num_runs
    [flux_0_obs(i,1),particle_counter_0(:,i),death_counter_0_obs(i,1),infection_counter_0_obs(i,1) ] = diffusion_with_obstacles_and_hit_count(2000, D, 0.01, 10.0, 0, 0, 100, penalty, bounding_box,3,new_particles);
end
fprintf('0 done');

flux_1_obs = zeros(num_runs,1);
particle_counter_1 = zeros(1000,num_runs); 
infection_counter_1_obs = zeros(num_runs,1);
death_counter_1_obs = zeros(num_runs,1);
for i = 1:num_runs
    [flux_1_obs(i,1),particle_counter_1(:,i),death_counter_1_obs(i,1),infection_counter_1_obs(i,1)] = diffusion_with_obstacles_and_hit_count(2000, D, 0.01, 10.0, 0, 1, 100, penalty, bounding_box,3,new_particles);
end
fprintf('1 done');

flux_5_obs = zeros(num_runs,1);
particle_counter_5 = zeros(1000,num_runs); 
infection_counter_5_obs = zeros(num_runs,1);
death_counter_5_obs = zeros(num_runs,1);
for i = 1:num_runs
    [flux_5_obs(i,1),particle_counter_5(:,i),death_counter_5_obs(i,1),infection_counter_5_obs(i,1)] = diffusion_with_obstacles_and_hit_count(2000, D, 0.01, 10.0, 0, 5, 100, penalty, bounding_box,3,new_particles);
end
fprintf('5 done');

flux_10_obs = zeros(num_runs,1);
particle_counter_10 = zeros(1000,num_runs); 
infection_counter_10_obs = zeros(num_runs,1);
death_counter_10_obs = zeros(num_runs,1);
for i = 1:num_runs
    [flux_10_obs(i,1),particle_counter_10(:,i),death_counter_10_obs(i,1),infection_counter_10_obs(i,1)] = diffusion_with_obstacles_and_hit_count(2000, D, 0.01, 10.0, 0, 10, 100, penalty, bounding_box,3,new_particles);
end
fprintf('10 done');

flux_50_obs = zeros(num_runs,1);
particle_counter_50 = zeros(1000,num_runs); 
infection_counter_50_obs = zeros(num_runs,1);
death_counter_50_obs = zeros(num_runs,1);
for i = 1:num_runs
    [flux_50_obs(i,1),particle_counter_50(:,i),death_counter_50_obs(i,1),infection_counter_50_obs(i,1)] = diffusion_with_obstacles_and_hit_count(2000, D, 0.01, 10.0, 0, 50, 100, penalty, bounding_box,3,new_particles);
end
fprintf('50 done');

flux_100_obs = zeros(num_runs,1);
particle_counter_100 = zeros(1000,num_runs); 
infection_counter_100_obs = zeros(num_runs,1);
death_counter_100_obs = zeros(num_runs,1);
for i = 1:num_runs
    [flux_100_obs(i,1),particle_counter_100(:,i),death_counter_100_obs(i,1),infection_counter_100_obs(i,1)] = diffusion_with_obstacles_and_hit_count(2000, D, 0.01, 10.0, 0, 100, 100, penalty, bounding_box,3,new_particles);
end
fprintf('100 done');


flux_500_obs = zeros(num_runs,1);
particle_counter_500 = zeros(1000,num_runs); 
infection_counter_500_obs = zeros(num_runs,1);
death_counter_500_obs = zeros(num_runs,1);
for i = 1:num_runs
    [flux_500_obs(i,1),particle_counter_500(:,i),death_counter_500_obs(i,1),infection_counter_500_obs(i,1)] = diffusion_with_obstacles_and_hit_count(2000, D, 0.01, 10.0, 0, 500, 100, penalty, bounding_box,3,new_particles);
end
fprintf('500 done');

flux_1000_obs = zeros(num_runs,1);
particle_counter_1000 = zeros(1000,num_runs); 
infection_counter_1000_obs = zeros(num_runs,1);
death_counter_1000_obs = zeros(num_runs,1);
for i = 1:num_runs
    [flux_1000_obs(i,1),particle_counter_1000(:,i),death_counter_1000_obs(i,1),infection_counter_1000_obs(i,1)] = diffusion_with_obstacles_and_hit_count(2000, D, 0.01, 10.0, 0, 1000, 100, penalty, bounding_box,3,new_particles);
end
fprintf('1000 done');
toc


figure;
boxplot([flux_0_obs,flux_1_obs,flux_5_obs,flux_10_obs,flux_50_obs,flux_100_obs,flux_500_obs,flux_1000_obs],...
    'Labels',{'0','1','5','10','50','100','500','1000'});
xlabel('Number of obstacles');
ylabel('Flux');


x=1:1000;
figure; 
plot_shaded_errorbar(x,particle_counter_0,'k','--'); hold on;
plot_shaded_errorbar(x,particle_counter_1,'r','--');
plot_shaded_errorbar(x,particle_counter_5,'g','--');
plot_shaded_errorbar(x,particle_counter_10,'b','--'); hold on;
plot_shaded_errorbar(x,particle_counter_50,'k'); hold on;
plot_shaded_errorbar(x,particle_counter_100,'r');
plot_shaded_errorbar(x,particle_counter_500,'g'); 
plot_shaded_errorbar(x,particle_counter_1000,'b'); hold on;
legend('0 obstacles','','1 obstacles','','5 obstacles','',...
    '10 obstacles','','50 obstacles','','100 obstacles','',...
    '500 obstacles','','1000 obstacles','');
ylabel('# particles');
xlabel('Time steps');

death_counter_vector = [mean(infection_counter_0_obs),...
    mean(death_counter_1_obs),mean(death_counter_5_obs),...
    mean(death_counter_10_obs),mean(death_counter_50_obs),...
    mean(death_counter_100_obs),mean(death_counter_500_obs),mean(death_counter_1000_obs)];
infection_counter_vector = [mean(infection_counter_0_obs),...
    mean(infection_counter_1_obs),mean(infection_counter_5_obs),...
    mean(infection_counter_10_obs),mean(infection_counter_50_obs),...
    mean(infection_counter_100_obs),mean(infection_counter_500_obs),mean(infection_counter_1000_obs)];
figure;
subplot(2,1,1)
scatter(1:8,infection_counter_vector);
xticklabels({'0','1','5','10','50','100','500','1000'});
xlabel('# obstacles');
ylabel('infection');
subplot(2,1,2)
scatter(1:8,death_counter_vector);
xticklabels({'0','1','5','10','50','100','500','1000'});
xlabel('# obstacles');
ylabel('Death');




obstacle_numbers = floor([2e3,1e3,1e3/(3^1),1e3/(3^2),1e3/(3^3),1e3/(3^4),1e3/(3^5),...
    1e3/(3^6),1e3/(3^7)]);
sialic_acid_penalty = [0:15:60];
mean_matrix = zeros(numel(obstacle_numbers),numel(sialic_acid_penalty));
std_matrix = zeros(numel(obstacle_numbers),numel(sialic_acid_penalty));
flux_matrix = zeros(10,numel(obstacle_numbers),numel(sialic_acid_penalty));
particle_counter_matrix = zeros(1000,10,numel(obstacle_numbers),numel(sialic_acid_penalty)); 
for num = 1:numel(obstacle_numbers) 
   for penalty = 1:numel(sialic_acid_penalty)
       for iteration = 1:10
           [flux_matrix(iteration,num,penalty),particle_counter_matrix(:,iteration,num,penalty)] = diffusion_with_obstacles_and_hit_count(2000, 1e6, 0.01, 10.0, 0, obstacle_numbers(num), 100, sialic_acid_penalty(penalty), [1e4,1e4,1e4],3,new_particles);
       end 
       fprintf('iteration over num= %d, penalty = %d\n',num,penalty);
       mean_matrix(num,penalty) = mean(flux_matrix(:,num,penalty));
       std_matrix(num,penalty) = std(flux_matrix(:,num,penalty));
   end
end

heatmap(obstacle_numbers,sialic_acid_penalty,mean_matrix');
% caxis([10,35]);
colormap turbo;
xlabel('# Obstacles');
ylabel('# Sialic Acids');

mean_matrix2 = mean_matrix;
hHM= heatmap(obstacle_numbers,sialic_acid_penalty,mean_matrix2'/max(mean_matrix2(:)));
% caxis([0.5,1]);
colormap turbo;
xlabel('# Obstacles');
ylabel('# Sialic Acids');
hHM.YDisplayData=flip(hHM.YDisplayData);  
hHM.XDisplayData=flip(hHM.XDisplayData);  
