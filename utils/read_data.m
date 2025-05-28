function [HIM,HIM_norm,GT,name] = read_data(data_ind)
    path='.\Dataset\';
    switch data_ind
        case 1
            %  Hydice urban
            inputs = 'HYDICE Urban cut';
            location = [path,inputs];
            load (location);
            inputs = 'HYDICE Urban cut gt';
            location = [path,inputs];
            load (location);
            name = {'HU'};
        case 2
            % Pavia
            inputs = 'pavia1011';
            location = [path,inputs];
            load (location);
            HIM = data;
            GT = double(map);
            name = {'Pavia'};
        case 3
            % Hyperion
            inputs = 'Hyperion_new';
            location = [path,inputs];
            load (location);
            HIM = Hyperion_new;
            inputs = 'Hyperion_gt';
            location = [path,inputs];
            load (location);
            GT = double(GT_Hyperion);
            name = {'Hyperion'};  

    end

        %% normalize for each band
        img = (ToVector(HIM))';
        [row,col,bands] = size(HIM);
        for i = 1:bands
            img(i,:) = (img(i,:) -min(img(i,:))) / (max(img(i,:))-min(img(i,:)));
        end
        HIM_norm = ToCube(img,row,col);

    
        %% figure show
        f_show = HIM(:,:,[16,43,72]);
        for i = 1:3
            max_f = max(max(f_show(:,:,i)));
            min_f = min(min(f_show(:,:,i)));
            f_show(:,:,i) = (f_show(:,:,i)-min_f)/(max_f-min_f);
        end
            figure(),imshow(f_show);

            figure(),imshow(double(GT),[]);
            num_GT = sum(sum(GT));
            fprintf('The number of GT is');  disp(num_GT); axis off



   