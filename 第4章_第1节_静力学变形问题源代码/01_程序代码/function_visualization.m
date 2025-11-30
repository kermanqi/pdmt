function[]=function_visualization(dimension,size_x,size_y,size_z,global_displacement_vector)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%一维情况
if dimension==1    
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
    figure(1)
    pcolor(displacement_x);shading interp;  % x方向位移结果可视化
    set(colorbar,'FontSize',14,'FontName','Times New Roman');
    colormap Jet
    
    displacement_y=zeros(size_y,size_x);    % y方向位移写为矩阵形式
    ii=0;                                   % 变量初始值取0
    for i=1:size_y                          % 按y方向粒子总数循环
        for j=1:size_x                      % 按x方向粒子总数循环
            ii=ii+1;                        % ii取值+1，从而遍历每个粒子
            displacement_y(i,j)=global_displacement_vector(ii*2,1);     % 总体位移向量中的y方向位移写入y方向位移矩阵
        end
    end
    figure(2)
    pcolor(displacement_y);shading interp;  % y方向位移结果可视化
    set(colorbar,'FontSize',14,'FontName','Times New Roman');
    colormap Jet
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%三维情况
if dimension==3
    [X,Y,Z]=meshgrid(1:size_z,1:size_y,1:size_x);
    displacement_x=zeros(size_z,size_y,size_x);    % x方向位移写为矩阵形式
    ii=0;                                   % 变量初始值取0
    for i=1:size_z                          % 按z方向粒子总数循环
        for j=1:size_y                      % 按y方向粒子总数循环
            for k=1:size_x                  % 按x方向粒子总数循环
                ii=ii+1;                        % ii取值+1，从而遍历每个粒子
                displacement_x(i,j,k)=global_displacement_vector(ii*3-2,1);   % 总体位移向量中的x方向位移写入x方向位移矩阵
            end
        end
    end
    figure(1)
    slice(X,Y,Z,displacement_x,[1,size_z],[1,size_y],[1,size_x]);  % x方向位移结果可视化
    set(gca,'DataAspectRatio',[1 1 1]);
    set(colorbar,'FontSize',14,'FontName','Times New Roman');
    colormap Jet
    
    displacement_y=zeros(size_z,size_y,size_x);    % y方向位移写为矩阵形式
    ii=0;                                   % 变量初始值取0
    for i=1:size_z                          % 按z方向粒子总数循环
        for j=1:size_y                      % 按y方向粒子总数循环
            for k=1:size_x                  % 按x方向粒子总数循环
                ii=ii+1;                        % ii取值+1，从而遍历每个粒子
                displacement_y(i,j,k)=global_displacement_vector(ii*3-1,1);   % 总体位移向量中的x方向位移写入y方向位移矩阵
            end
        end
    end
    figure(2)
    slice(X,Y,Z,displacement_y,[1,size_z],[1,size_y],[1,size_x]);  % y方向位移结果可视化
    set(gca,'DataAspectRatio',[1 1 1]);
    set(colorbar,'FontSize',14,'FontName','Times New Roman');
    colormap Jet
    
    displacement_z=zeros(size_z,size_y,size_x);    % z方向位移写为矩阵形式
    ii=0;                                   % 变量初始值取0
    for i=1:size_z                          % 按z方向粒子总数循环
        for j=1:size_y                      % 按y方向粒子总数循环
            for k=1:size_x                  % 按x方向粒子总数循环
                ii=ii+1;                        % ii取值+1，从而遍历每个粒子
                displacement_z(i,j,k)=global_displacement_vector(ii*3,1);   % 总体位移向量中的x方向位移写入z方向位移矩阵
            end
        end
        
        
        
    end
    figure(3)
    slice(X,Y,Z,displacement_z,[1,size_z],[1,size_y],[1,size_x]);  % z方向位移结果可视化
    set(gca,'DataAspectRatio',[1 1 1]);
    set(colorbar,'FontSize',14,'FontName','Times New Roman');
    colormap Jet

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%