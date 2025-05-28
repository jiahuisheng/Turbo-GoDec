function a = plot_image(detection_cub,GT,data_name,method)

[N,w] = size(detection_cub); 
[row,col] = size(GT);


for i = 1:w 
    detection = reshape(detection_cub(:,i),[row,col]); 
    method_name = method{i};
    figure(),colormap;imagesc(detection);axis off;%axis equal

    % 自动调整绘图窗口大小并去除多余的边距
    set(gcf, 'Units', 'pixels');
    position = get(gcf, 'Position');
    aspect_ratio = size(detection, 2) / size(detection, 1);
    new_width = position(4) * aspect_ratio;
    set(gcf, 'Position', [position(1), position(2), new_width, position(4)]);
    % set (gcf,'Position',[0,0,512,512]);
    axis image; set(gca,'xtick',[]); set(gca,'ytick',[]);set(gca,'Position',[0 0 1 1]);

    % frame = getframe;
    % color_image = frame2im(frame);
    folder = 'figures';
    extension = '.jpg';
    % imwrite(color_image, fullfile(folder, [char(data_name),'_',method_name,extension]));
    saveas(gcf, fullfile(folder, [char(data_name),'_',method_name,extension]), 'jpg');
    close(gcf);

end
a = w;