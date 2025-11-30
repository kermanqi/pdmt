function[]=function_visualization(dimension,size_x,size_y,size_z,global_displacement_vector,boundary_condition_left,boundary_condition_right, ...
boundary_condition_bottom,boundary_condition_top,boundary_condition_behind,boundary_condition_front,m)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%一维情况
if dimension==1   
    if boundary_condition_right(1)==1             % 如果右边界施加位移条件，则执行条件语句中的内容
        global_displacement_vector(size_x-m+1:size_x)=[];
    end
    if boundary_condition_left(1)==1              % 如果左边界施加位移条件，则执行条件语句中的内容
        global_displacement_vector(1:m)=[];
    end
    figure(1)
    plot(global_displacement_vector);  % x方向位移结果可视化
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%二维情况
if dimension==2
    displacement_x=zeros(size_y,size_x);    % x方向位移写为矩阵形式
    ii=0;                                   % 变量初始值取0
    for i=1:size_y                          % 按y方向粒子总数循环
        for j=1:size_x                      % 按x方向粒子总数循环
            ii=ii+1;                        % ii取值+1，从而遍历每个粒子
            displacement_x(i,j)=global_displacement_vector(ii*2-1,1);   % 总体位移向量中的x方向位移写入x方向位移矩阵
        end
    end
    displacement_y=zeros(size_y,size_x);    % y方向位移写为矩阵形式
    ii=0;                                   % 变量初始值取0
    for i=1:size_y                          % 按y方向粒子总数循环
        for j=1:size_x                      % 按x方向粒子总数循环
            ii=ii+1;                        % ii取值+1，从而遍历每个粒子
            displacement_y(i,j)=global_displacement_vector(ii*2,1);     % 总体位移向量中的y方向位移写入y方向位移矩阵
        end
    end
    
    if boundary_condition_right(1)==1 || boundary_condition_right(2)==1               % 如果右边界施加位移条件，则执行条件语句中的内容
        displacement_x(:,size_x-m+1:size_x)=[];
        displacement_y(:,size_x-m+1:size_x)=[];
    end    
    if boundary_condition_left(1)==1 || boundary_condition_left(2)==1               % 如果左边界施加位移条件，则执行条件语句中的内容
        displacement_x(:,1:m)=[];
        displacement_y(:,1:m)=[];
    end
    if boundary_condition_top(1)==1 || boundary_condition_top(2)==1            % 如果上边界施加位移条件，则执行条件语句中的内容
        displacement_x(size_y-m+1:size_y,:)=[];
        displacement_y(size_y-m+1:size_y,:)=[];
    end
    if boundary_condition_bottom(1)==1 || boundary_condition_bottom(2)==1         % 如果下边界施加位移条件，则执行条件语句中的内容
        displacement_x(1:m,:)=[];
        displacement_y(1:m,:)=[];
    end

    figure(1)
    pcolor(displacement_x);shading interp;  % x方向位移结果可视化
    set(colorbar,'FontSize',14,'FontName','Times New Roman');  % 设置图例的字号字体
    colormap Jet  % 图像颜色映射设置为Jet
    
    figure(2)
    pcolor(displacement_y);shading interp;  % y方向位移结果可视化
    set(colorbar,'FontSize',14,'FontName','Times New Roman');  % 设置图例的字号字体
    colormap Jet  % 图像颜色映射设置为Jet
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%三维情况
if dimension==3
    displacement_x=zeros(size_x,size_y,size_z);    % x方向位移写为矩阵形式
    ii=0;                                   % 变量初始值取0
    for i=1:size_z                          % 按z方向粒子总数循环
        for j=1:size_y                      % 按y方向粒子总数循环
            for k=1:size_x                  % 按x方向粒子总数循环
                ii=ii+1;                        % ii取值+1，从而遍历每个粒子
                displacement_x(k,j,i)=global_displacement_vector(ii*3-2,1);   % 总体位移向量中的x方向位移写入x方向位移矩阵
            end
        end
    end
    displacement_y=zeros(size_x,size_y,size_z);    % y方向位移写为矩阵形式
    ii=0;                                   % 变量初始值取0
    for i=1:size_z                          % 按z方向粒子总数循环
        for j=1:size_y                      % 按y方向粒子总数循环
            for k=1:size_x                  % 按x方向粒子总数循环
                ii=ii+1;                        % ii取值+1，从而遍历每个粒子
                displacement_y(k,j,i)=global_displacement_vector(ii*3-1,1);   % 总体位移向量中的y方向位移写入y方向位移矩阵
            end
        end
    end
    displacement_z=zeros(size_x,size_y,size_z);    % z方向位移写为矩阵形式
    ii=0;                                   % 变量初始值取0
    for i=1:size_z                          % 按z方向粒子总数循环
        for j=1:size_y                      % 按y方向粒子总数循环
            for k=1:size_x                  % 按x方向粒子总数循环
                ii=ii+1;                        % ii取值+1，从而遍历每个粒子
                displacement_z(k,j,i)=global_displacement_vector(ii*3,1);   % 总体位移向量中的z方向位移写入z方向位移矩阵
            end
        end
    end
    
    size_boundary_x=0;                      % 用来统计x方向施加位移边界条件的个数，该变量初始值取0
    size_boundary_y=0;                      % 用来统计y方向施加位移边界条件的个数，该变量初始值取0
    size_boundary_z=0;                      % 用来统计z方向施加位移边界条件的个数，该变量初始值取0   
    if boundary_condition_right(1)==1 || boundary_condition_right(2)==1 || boundary_condition_right(3)==1               % 如果右边界施加位移条件，则执行条件语句中的内容
        size_boundary_x=size_boundary_x+1;      % x方向施加位移边界条件的个数+1
        displacement_x(size_x-m+1:size_x,:,:)=[];
        displacement_y(size_x-m+1:size_x,:,:)=[];
        displacement_z(size_x-m+1:size_x,:,:)=[];
    end    
    if boundary_condition_left(1)==1 || boundary_condition_left(2)==1 || boundary_condition_left(3)==1                % 如果左边界施加位移条件，则执行条件语句中的内容
        size_boundary_x=size_boundary_x+1;      % x方向施加位移边界条件的个数+1
        displacement_x(1:m,:,:)=[];
        displacement_y(1:m,:,:)=[];
        displacement_z(1:m,:,:)=[];
    end
    if boundary_condition_top(1)==1 || boundary_condition_top(2)==1 || boundary_condition_top(3)==1            % 如果上边界施加位移条件，则执行条件语句中的内容
        size_boundary_y=size_boundary_y+1;  % y方向施加位移边界条件的个数+1
        displacement_x(:,size_y-m+1:size_y,:)=[];
        displacement_y(:,size_y-m+1:size_y,:)=[];
        displacement_z(:,size_y-m+1:size_y,:)=[];
    end
    if boundary_condition_bottom(1)==1 || boundary_condition_bottom(2)==1 || boundary_condition_bottom(3)==1         % 如果下边界施加位移条件，则执行条件语句中的内容
        size_boundary_y=size_boundary_y+1;  % y方向施加位移边界条件的个数+1
        displacement_x(:,1:m,:)=[];
        displacement_y(:,1:m,:)=[];
        displacement_z(:,1:m,:)=[];
    end
    if boundary_condition_front(1)==1 || boundary_condition_front(2)==1 || boundary_condition_front(3)==1            % 如果前边界施加位移条件，则执行条件语句中的内容
        size_boundary_z=size_boundary_z+1;  % y方向施加位移边界条件的个数+1
        displacement_x(:,:,size_z-m+1:size_z)=[];
        displacement_y(:,:,size_z-m+1:size_z)=[];
        displacement_z(:,:,size_z-m+1:size_z)=[];
    end
    if boundary_condition_behind(1)==1 || boundary_condition_behind(2)==1 || boundary_condition_behind(3)==1         % 如果后边界施加位移条件，则执行条件语句中的内容
        size_boundary_z=size_boundary_z+1;  % y方向施加位移边界条件的个数+1
        displacement_x(:,:,1:m)=[];
        displacement_y(:,:,1:m)=[];
        displacement_z(:,:,1:m)=[];
    end

    size_x=size_x-size_boundary_x*m;
    size_y=size_y-size_boundary_y*m;
    size_z=size_z-size_boundary_z*m;
    
    figure(1)
    slice(displacement_x,[1,size_y],[1,size_x],[1,size_z]);  % x方向位移结果可视化
    set(gca,'DataAspectRatio',[1 1 1]);  % 设置坐标轴比例
    axis off
    set(colorbar,'FontSize',14,'FontName','Times New Roman');  % 设置图例的字号字体
    colormap Jet  % 图像颜色映射设置为Jet
    
    figure(2)
    slice(displacement_y,[1,size_y],[1,size_x],[1,size_z]);  % y方向位移结果可视化
    set(gca,'DataAspectRatio',[1 1 1]);  % 设置坐标轴比例
    axis off
    set(colorbar,'FontSize',14,'FontName','Times New Roman');  % 设置图例的字号字体
    colormap Jet  % 图像颜色映射设置为Jet
    
    figure(3)
    slice(displacement_z,[1,size_y],[1,size_x],[1,size_z]);  % z方向位移结果可视化
    set(gca,'DataAspectRatio',[1 1 1]);  % 设置坐标轴比例
    axis off
    set(colorbar,'FontSize',14,'FontName','Times New Roman');  % 设置图例的字号字体
    colormap Jet  % 图像颜色映射设置为Jet

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%