function [outputArg1,outputArg2] = plot_shaded_errorbar(x,y,c,sign)




mean_y = mean(y,2);
ste_y = std(y,0,2)/sqrt(size(y,2));
if nargin==3
    plot(x,mean_y,'Color',c,'LineWidth',2); hold on;
else
    plot(x,mean_y,sign,'Color',c,'LineWidth',2); hold on;
end

patch([x'; flip(x')], [mean_y-ste_y; flip(mean_y+ste_y)],...
    c,'FaceAlpha',0.1,'EdgeColor','none');
end

