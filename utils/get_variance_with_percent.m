function variance_percent = get_variance_with_percent(x,percent,sort_flag)
    % 1. 计算向量 x 的绝对值
    abs_x = abs(x);
    
    % 2. 找出绝对值中前百分之10的元素
    n = length(x);
    percent_idx = round(percent * n);
    sorted_abs_x = sort(abs_x, sort_flag);
    percent_values = sorted_abs_x(1:percent_idx);
    
    % 3. 计算这些元素的方差
    variance_percent = var(percent_values);
end