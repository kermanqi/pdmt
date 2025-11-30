function [particle,size_particle,size_x,size_y,size_z,size_displacement_boundary_left,particle_boundary_left,size_displacement_boundary_right,particle_boundary_right, ...
          size_displacement_boundary_bottom,particle_boundary_bottom,size_displacement_boundary_top,particle_boundary_top,size_displacement_boundary_behind, ...
          particle_boundary_behind,size_displacement_boundary_front,particle_boundary_front,size_uniform_load_left,particle_uniform_load_left,size_uniform_load_right, ...
          particle_uniform_load_right,size_uniform_load_bottom,particle_uniform_load_bottom,size_uniform_load_top,particle_uniform_load_top,size_uniform_load_behind, ...
          particle_uniform_load_behind,size_uniform_load_front,particle_uniform_load_front, ...
          size_temperature_boundary_left,particle_boundary_left_temperature,size_temperature_boundary_right,particle_boundary_right_temperature, ...
          size_temperature_boundary_bottom,particle_boundary_bottom_temperature,size_temperature_boundary_top,particle_boundary_top_temperature,size_temperature_boundary_behind, ...
          particle_boundary_behind_temperature,size_temperature_boundary_front,particle_boundary_front_temperature,size_heat_flux_left,particle_heat_flux_left,size_heat_flux_right, ...
          particle_heat_flux_right,size_heat_flux_bottom,particle_heat_flux_bottom,size_heat_flux_top,particle_heat_flux_top,size_heat_flux_behind, ...
          particle_heat_flux_behind,size_heat_flux_front,particle_heat_flux_front,global_velocity_vector,global_displacement_vector,global_temperature_vector]=function_discretization(dimension,boundary_condition_left, ...
          boundary_condition_right,boundary_condition_bottom,boundary_condition_top,boundary_condition_behind,boundary_condition_front,xx_min,yy_min,zz_min,xx_max,yy_max,zz_max,length_x,length_y,length_z,dettx,m, ...
          boundary_condition_left_temperature,boundary_condition_right_temperature, ...
    boundary_condition_bottom_temperature,boundary_condition_top_temperature,boundary_condition_behind_temperature,boundary_condition_front_temperature,temperature_initial)

size_displacement_boundary_left=0;  % 左边界施加位移边界条件的粒子总数赋初值
particle_boundary_left=[];          % 定义左边界施加位移边界条件的粒子集合为空集
size_displacement_boundary_right=0; % 右边界施加位移边界条件的粒子总数赋初值
particle_boundary_right=[];         % 定义右边界施加位移边界条件的粒子集合为空集
size_displacement_boundary_bottom=0; % 下边界施加位移边界条件的粒子总数赋初值
particle_boundary_bottom=[];        % 定义下边界施加位移边界条件的粒子集合为空集
size_displacement_boundary_top=0;   % 上边界施加位移边界条件的粒子总数赋初值
particle_boundary_top=[];           % 定义上边界施加位移边界条件的粒子集合为空集
size_displacement_boundary_behind=0; % 后边界施加位移边界条件的粒子总数赋初值
particle_boundary_behind=[];           % 定义后边界施加位移边界条件的粒子集合为空集
size_displacement_boundary_front=0; % 前边界施加位移边界条件的粒子总数赋初值
particle_boundary_front=[];         % 定义前边界施加位移边界条件的粒子集合为空集

size_uniform_load_left=0;       % 左边界施加均布荷载条件的粒子总数赋初值
particle_uniform_load_left=[];  % 定义左边界施加均布荷载条件的粒子集合为空集
size_uniform_load_right=0;      % 右边界施加均布荷载条件的粒子总数赋初值
particle_uniform_load_right=[]; % 定义右边界施加均布荷载条件的粒子集合为空集
size_uniform_load_bottom=0;     % 下边界施加均布荷载条件的粒子总数赋初值
particle_uniform_load_bottom=[];% 定义下边界施加均布荷载条件的粒子集合为空集
size_uniform_load_top=0;        % 上边界施加均布荷载条件的粒子总数赋初值
particle_uniform_load_top=[];   % 定义上边界施加均布荷载条件的粒子集合为空集
size_uniform_load_behind=0;     % 后边界施加均布荷载条件的粒子总数赋初值
particle_uniform_load_behind=[];% 定义后边界施加均布荷载条件的粒子集合为空集
size_uniform_load_front=0;      % 前边界施加均布荷载条件的粒子总数赋初值
particle_uniform_load_front=[]; % 定义前边界施加均布荷载条件的粒子集合为空集
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%一维情况
if dimension==1    
    size_boundary_x=0;                      % 用来统计x方向施加位移边界条件的个数，该变量初始值取0
    size_boundary_left=0;                   % 用于表示左边界是否施加位移条件，=0表示左边界未施加位移条件，=1表示左边界施加位移条件，该变量初始值取0
    if boundary_condition_left(1)==1 || boundary_condition_left_temperature==1             % 如果左边界施加位移条件，则执行条件语句中的内容
        size_boundary_x=size_boundary_x+1;      % x方向施加位移边界条件的个数+1
        size_boundary_left=1;                   % =1 代表左边界施加位移条件
    end
    if boundary_condition_right(1)==1 || boundary_condition_right_temperature==1             % 如果右边界施加位移条件，则执行条件语句中的内容
        size_boundary_x=size_boundary_x+1;      % x方向施加位移边界条件的个数+1
    end
    
    size_x=length_x/dettx+size_boundary_x*m;     % x方向粒子总数=模型x方向长度/每个粒子长度+x方向施加位移边界条件的个数*粒子作用的层数
    size_y=1;     % 对于一维问题，y方向粒子总数=1
    size_z=1;     % 对于一维问题，z方向粒子总数=1
    ii=0;       % 变量初始值取0
    for i=1:size_x          % 按x方向粒子总数循环
        ii=ii+1;        % ii取值+1，从而遍历每个粒子
        particle(ii,1)=xx_min+(i-0.5)*dettx-size_boundary_left*m*dettx;     % 第ii个粒子的x方向坐标=x方向最小坐标*(i-0.5)*每个粒子长度，如果左边界定义了位移边界，第ii个粒子的x方向坐标=x方向最小坐标*(i-0.5)*每个粒子长度-粒子作用的层数*每个粒子长度
    end
    size_particle=size(particle,1);          % 粒子总数
%%%%%%%%%%%%%%%%%%%%%
% 施加位移边界条件的粒子
%左边界
    ii_displacement_boundary_left=0;    % 左边界施加位移边界条件的粒子计数赋初值为0
    if boundary_condition_left(1)==1    % 如果左边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            if xx_i<xx_min              % 如果xx_i小于x方向最小坐标，说明粒子i处于左边界，则执行条件语句中的内容
                ii_displacement_boundary_left=ii_displacement_boundary_left+1;  % 左边界施加位移边界条件粒子计数加1
                particle_boundary_left(ii_displacement_boundary_left)=i;        % 在集合particle_boundary_left中记录粒子i的标号
            end
        end
        size_displacement_boundary_left=ii_displacement_boundary_left;  % 记录左边界施加位移边界条件的粒子总数
    end
%右边界
    ii_displacement_boundary_right=0;    % 右边界施加位移边界条件的粒子计数赋初值为0
    if boundary_condition_right(1)==1    % 如果右边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            if xx_i>xx_max              % 如果xx_i大于x方向最大坐标，说明粒子i处于右边界，则执行条件语句中的内容
                ii_displacement_boundary_right=ii_displacement_boundary_right+1;  % 右边界施加位移边界条件粒子计数加1
                particle_boundary_right(ii_displacement_boundary_right)=i;        % 在集合particle_boundary_right中记录粒子i的标号
            end
        end
        size_displacement_boundary_right=ii_displacement_boundary_right;  % 记录右边界施加位移边界条件的粒子总数
    end
%%%%%%%%%%%%%%%%%%%%%
% 施加均布荷载条件的粒子
%左边界
    ii_uniform_load_left=0;                                         % 左边界施加均布荷载条件的粒子计数赋初值为0
    if boundary_condition_left(1)==2                                % 如果左边界施加均布荷载条件，则执行条件语句中的内容
        for i=1:size_particle                                       % 对所有粒子进行循环
            xx_i=particle(i,1);                                     % xx_i表示粒子i在x方向的坐标
            if xx_i<(xx_min+dettx) && xx_i>xx_min                                  % 如果xx_i小于(x方向最小坐标+dettx)，说明粒子i为位于最左侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_left=ii_uniform_load_left+1;        % 左边界施加均布荷载条件粒子计数加1
                particle_uniform_load_left(ii_uniform_load_left)=i; % 在集合particle_uniform_load_left中记录粒子i的标号
            end
        end
        size_uniform_load_left=ii_uniform_load_left;                % 记录左边界施加均布荷载条件的粒子总数
    end
%右边界
    ii_uniform_load_right=0;                                        % 右边界施加均布荷载条件的粒子计数赋初值为0
    if boundary_condition_right(1)==2                               % 如果右边界施加均布荷载条件，则执行条件语句中的内容 
        for i=1:size_particle                                       % 对所有粒子进行循环
            xx_i=particle(i,1);                                     % xx_i表示粒子i在x方向的坐标
            if xx_i>(xx_max-dettx) && xx_i<xx_max                                  % 如果xx_i大于(x方向最大坐标-dettx)，说明粒子i为位于最右侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_right=ii_uniform_load_right+1;        % 右边界施加均布荷载条件粒子计数加1
                particle_uniform_load_right(ii_uniform_load_right)=i; % 在集合particle_uniform_load_right中记录粒子i的标号
        	end
        end
        size_uniform_load_right=ii_uniform_load_right;                % 记录右边界施加均布荷载条件的粒子总数
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%二维情况
if dimension==2
    size_boundary_x=0;                      % 用来统计x方向施加位移边界条件的个数，该变量初始值取0
    size_boundary_left=0;                   % 用于表示左边界是否施加位移条件，=0表示左边界未施加位移条件，=1表示左边界施加位移条件，该变量初始值取0
    if boundary_condition_left(1)==1 || boundary_condition_left(2)==1 || boundary_condition_left_temperature==1               % 如果左边界施加位移条件，则执行条件语句中的内容
        size_boundary_x=size_boundary_x+1;      % x方向施加位移边界条件的个数+1
        size_boundary_left=1;                   % =1 代表左边界施加位移条件
    end
    if boundary_condition_right(1)==1 || boundary_condition_right(2)==1 || boundary_condition_right_temperature==1               % 如果右边界施加位移条件，则执行条件语句中的内容
        size_boundary_x=size_boundary_x+1;      % x方向施加位移边界条件的个数+1
    end
    size_boundary_y=0;                      % 用来统计y方向施加位移边界条件的个数，该变量初始值取0
    size_boundary_bottom=0;                 % 用于表示下边界是否施加位移条件，=0表示下边界未施加位移条件，=1表示下边界施加位移条件，该变量初始值取0
    if boundary_condition_bottom(1)==1 || boundary_condition_bottom(2)==1 || boundary_condition_bottom_temperature==1         % 如果下边界施加位移条件，则执行条件语句中的内容
        size_boundary_y=size_boundary_y+1;  % y方向施加位移边界条件的个数+1
        size_boundary_bottom=1;              % =1 代表下边界施加位移条件
    end
    if boundary_condition_top(1)==1 || boundary_condition_top(2)==1 || boundary_condition_top_temperature==1            % 如果上边界施加位移条件，则执行条件语句中的内容
        size_boundary_y=size_boundary_y+1;  % y方向施加位移边界条件的个数+1
    end

    size_x=length_x/dettx+size_boundary_x*m;     % x方向粒子总数=模型x方向长度/每个粒子长度+x方向施加位移边界条件的个数*粒子作用的层数
    size_y=length_y/dettx+size_boundary_y*m;     % y方向粒子总数=模型y方向长度/每个粒子长度+y方向施加位移边界条件的个数*粒子作用的层数
    size_z=1;                                           % 对于二维问题，z方向粒子总数=1
    ii=0;       % 变量初始值取0
    for i=1:size_y          % 按y方向粒子总数循环
        for j=1:size_x      % 按x方向粒子总数循环
            ii=ii+1;        % ii取值+1，从而遍历每个粒子
            particle(ii,1)=xx_min+(j-0.5)*dettx-size_boundary_left*m*dettx;     % 第ii个粒子的x方向坐标=x方向最小坐标*(j-0.5)*每个粒子长度，如果左边界定义了位移边界，第ii个粒子的x方向坐标=x方向最小坐标*(j-0.5)*每个粒子长度-粒子作用的层数*每个粒子长度
            particle(ii,2)=yy_min+(i-0.5)*dettx-size_boundary_bottom*m*dettx;	% 第ii个粒子的y方向坐标=y方向最小坐标*(i-0.5)*每个粒子长度，如果下边界定义了位移边界，第ii个粒子的y方向坐标=y方向最小坐标*(i-0.5)*每个粒子长度-粒子作用的层数*每个粒子长度
        end
    end
    size_particle=size(particle,1);          % 粒子总数

%%%%%%%%%%%%%%%%%%%%%
% 施加位移边界条件的粒子
%左边界
    ii_displacement_boundary_left=0; % 左边界施加位移边界条件粒子计数赋初值
    if boundary_condition_left(1)==1 || boundary_condition_left(2)==1 % 如果左边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle                                         % 对所有粒子进行循环
            xx_i=particle(i,1);                                       % xx_i表示粒子i在x方向的坐标
            if xx_i<xx_min                                            % 如果xx_i小于x方向最小坐标，说明粒子i处于左边界，则执行条件语句中的内容
                ii_displacement_boundary_left=ii_displacement_boundary_left+1;    % 左边界施加位移边界条件粒子计数加1
                particle_boundary_left(ii_displacement_boundary_left)=i;          % 在集合particle_boundary_left中记录粒子i的标号
            end
        end
        size_displacement_boundary_left=ii_displacement_boundary_left;            % 记录左边界施加位移边界条件的粒子总数
    end
%右边界
    ii_displacement_boundary_right=0;   % 右边界施加位移边界条件粒子计数赋初值
    if boundary_condition_right(1)==1 || boundary_condition_right(2)==1 % 如果右边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle                                           % 对所有粒子进行循环
            xx_i=particle(i,1);                                         % xx_i表示粒子i在x方向的坐标
            if xx_i>xx_max                                              % 如果xx_i大于x方向最大坐标，说明粒子i处于右边界，则执行条件语句中的内容
                ii_displacement_boundary_right=ii_displacement_boundary_right+1;   % 右边界施加位移边界条件粒子计数加1
                particle_boundary_right(ii_displacement_boundary_right)=i;         % 在集合particle_boundary_right中记录粒子i的标号
            end
        end
        size_displacement_boundary_right=ii_displacement_boundary_right;           % 记录右边界施加位移边界条件的粒子总数
    end
%下边界
    ii_displacement_boundary_bottom=0;   % 下边界施加位移边界条件粒子计数赋初值
    if boundary_condition_bottom(1)==1 || boundary_condition_bottom(2)==1 % 如果下边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle                                             % 对所有粒子进行循环
            yy_i=particle(i,2);                                           % yy_i表示粒子i在y方向的坐标
            if yy_i<yy_min                                                % 如果yy_i小于y方向最小坐标，说明粒子i处于下边界，则执行条件语句中的内容
                ii_displacement_boundary_bottom=ii_displacement_boundary_bottom+1; % 下边界施加位移边界条件粒子计数加1
                particle_boundary_bottom(ii_displacement_boundary_bottom)=i;       % 在集合particle_boundary_bottom中记录粒子i的标号
            end
        end
        size_displacement_boundary_bottom=ii_displacement_boundary_bottom;         % 记录下边界施加位移边界条件的粒子总数
    end
%上边界
    ii_displacement_boundary_top=0;   % 上边界施加位移边界条件粒子计数赋初值
    if boundary_condition_top(1)==1 || boundary_condition_top(2)==1      % 如果上边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle                                            % 对所有粒子进行循环
            yy_i=particle(i,2);                                          % yy_i表示粒子i在y方向的坐标
            if yy_i>yy_max                                               % 如果yy_i大于y方向最大坐标，说明粒子i处于上边界，则执行条件语句中的内容
                ii_displacement_boundary_top=ii_displacement_boundary_top+1;       % 上边界施加位移边界条件粒子计数加1
                particle_boundary_top(ii_displacement_boundary_top)=i;             % 在集合particle_boundary_top中记录粒子i的标号
            end
        end
        size_displacement_boundary_top=ii_displacement_boundary_top;               % 记录上边界施加位移边界条件的粒子总数
    end
%%%%%%%%%%%%%%%%%%%%%
% 施加均布荷载条件的粒子
%左边界
    ii_uniform_load_left=0; % 左边界施加位移边界条件粒子计数赋初值
    if boundary_condition_left(1)==2 || boundary_condition_left(2)==2     % 如果左边界施加均布荷载条件，则执行条件语句中的内容
        for i=1:size_particle                                             % 对所有粒子进行循环
            xx_i=particle(i,1);                                           % xx_i表示粒子i在x方向的坐标
            if xx_i<(xx_min+dettx) && xx_i>xx_min                                        % 如果xx_i小于(x方向最小坐标+dettx)，说明粒子i为位于最左侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_left=ii_uniform_load_left+1;              % 左边界施加均布荷载条件粒子计数加1
                particle_uniform_load_left(ii_uniform_load_left)=i;       % 在集合particle_uniform_load_left中记录粒子i的标号
            end
        end
        size_uniform_load_left=ii_uniform_load_left;                      % 记录左边界施加均布荷载条件的粒子总数
    end
%右边界
    ii_uniform_load_right=0;   % 右边界施加位移边界条件粒子计数赋初值
    if boundary_condition_right(1)==2 || boundary_condition_right(2)==2   % 如果右边界施加均布荷载条件，则执行条件语句中的内容 
        for i=1:size_particle                                             % 对所有粒子进行循环
            xx_i=particle(i,1);                                           % xx_i表示粒子i在x方向的坐标
            if xx_i>(xx_max-dettx) && xx_i<xx_max                                        % 如果xx_i大于(x方向最大坐标-dettx)，说明粒子i为位于最右侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_right=ii_uniform_load_right+1;            % 右边界施加均布荷载条件粒子计数加1
                particle_uniform_load_right(ii_uniform_load_right)=i;     % 在集合particle_uniform_load_right中记录粒子i的标号
            end
        end
        size_uniform_load_right=ii_uniform_load_right;                    % 记录右边界施加均布荷载条件的粒子总数
    end
%下边界
    ii_uniform_load_bottom=0;   % 下边界施加位移边界条件粒子计数赋初值
    if boundary_condition_bottom(1)==2 || boundary_condition_bottom(2)==2 % 如果下边界施加均布荷载条件，则执行条件语句中的内容
        for i=1:size_particle                                             % 对所有粒子进行循环
            yy_i=particle(i,2);                                           % yy_i表示粒子i在y方向的坐标
            if yy_i<(yy_min+dettx) && yy_i>yy_min                                        % 如果yy_i小于(y方向最小坐标+dettx)，说明粒子i为位于最下侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_bottom=ii_uniform_load_bottom+1;          % 下边界施加均布荷载条件粒子计数加1
                particle_uniform_load_bottom(ii_displacement_boundary_bottom)=i;    % 在集合particle_uniform_load_bottom中记录粒子i的标号
            end
        end
        size_uniform_load_bottom=ii_displacement_uniform_load_bottom;     % 记录下边界施加均布荷载条件的粒子总数
    end
%上边界
    ii_uniform_load_top=0;   % 上边界施加位移边界条件粒子计数赋初值
    if boundary_condition_top(1)==2 || boundary_condition_top(2)==2       % 如果上边界施加均布荷载条件，则执行条件语句中的内容
        for i=1:size_particle                                             % 对所有粒子进行循环
            yy_i=particle(i,2);                                           % yy_i表示粒子i在y方向的坐标
            if yy_i>(yy_max-dettx) && yy_i<yy_max                                        % 如果yy_i大于(y方向最大坐标-dettx)，说明粒子i为位于最上侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_top=ii_uniform_load_top+1;                % 上边界施加均布荷载条件粒子计数加1
                particle_uniform_load_top(ii_uniform_load_top)=i;         % 在集合particle_uniform_load_top中记录粒子i的标号
            end
        end
        size_uniform_load_top=ii_uniform_load_top;                        % 记录上边界施加均布荷载条件的粒子总数
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%三维情况
if dimension==3
    size_boundary_x=0;                      % 用来统计x方向施加位移边界条件的个数，该变量初始值取0
    size_boundary_left=0;                   % 用于表示左边界是否施加位移条件，=0表示左边界未施加位移条件，=1表示左边界施加位移条件，该变量初始值取0
    if boundary_condition_left(1)==1 || boundary_condition_left(2)==1 || boundary_condition_left(3)==1 || boundary_condition_left_temperature==1                % 如果左边界施加位移条件，则执行条件语句中的内容
        size_boundary_x=size_boundary_x+1;      % x方向施加位移边界条件的个数+1
        size_boundary_left=1;                   % =1 代表左边界施加位移条件
    end
    if boundary_condition_right(1)==1 || boundary_condition_right(2)==1 || boundary_condition_right(3)==1 || boundary_condition_right_temperature==1               % 如果右边界施加位移条件，则执行条件语句中的内容
        size_boundary_x=size_boundary_x+1;      % x方向施加位移边界条件的个数+1
    end
    size_boundary_y=0;                      % 用来统计y方向施加位移边界条件的个数，该变量初始值取0
    size_boundary_bottom=0;                 % 用于表示下边界是否施加位移条件，=0表示下边界未施加位移条件，=1表示下边界施加位移条件，该变量初始值取0
    if boundary_condition_bottom(1)==1 || boundary_condition_bottom(2)==1 || boundary_condition_bottom(3)==1 || boundary_condition_bottom_temperature==1         % 如果下边界施加位移条件，则执行条件语句中的内容
        size_boundary_y=size_boundary_y+1;  % y方向施加位移边界条件的个数+1
        size_boundary_bottom=1;              % =1 代表下边界施加位移条件
    end
    if boundary_condition_top(1)==1 || boundary_condition_top(2)==1 || boundary_condition_top(3)==1 || boundary_condition_top_temperature==1            % 如果上边界施加位移条件，则执行条件语句中的内容
        size_boundary_y=size_boundary_y+1;  % y方向施加位移边界条件的个数+1
    end
    size_boundary_z=0;                      % 用来统计z方向施加位移边界条件的个数，该变量初始值取0
    size_boundary_behind=0;                 % 用于表示后边界是否施加位移条件，=0表示后边界未施加位移条件，=1表示后边界施加位移条件，该变量初始值取0
    if boundary_condition_behind(1)==1 || boundary_condition_behind(2)==1 || boundary_condition_behind(3)==1 || boundary_condition_behind_temperature==1         % 如果后边界施加位移条件，则执行条件语句中的内容
        size_boundary_z=size_boundary_z+1;  % y方向施加位移边界条件的个数+1
        size_boundary_behind=1;              % =1 代表下边界施加位移条件
    end
    if boundary_condition_front(1)==1 || boundary_condition_front(2)==1 || boundary_condition_front(3)==1 || boundary_condition_front_temperature==1            % 如果前边界施加位移条件，则执行条件语句中的内容
        size_boundary_z=size_boundary_z+1;  % y方向施加位移边界条件的个数+1
    end
    size_x=length_x/dettx+size_boundary_x*m;     % x方向粒子总数=模型x方向长度/每个粒子长度+x方向施加位移边界条件的个数*粒子作用的层数
    size_y=length_y/dettx+size_boundary_y*m;     % y方向粒子总数=模型y方向长度/每个粒子长度+y方向施加位移边界条件的个数*粒子作用的层数
    size_z=length_z/dettx+size_boundary_z*m;     % z方向粒子总数=模型y方向长度/每个粒子长度+y方向施加位移边界条件的个数*粒子作用的层数
    ii=0;       % 变量初始值取0
    for i=1:size_z              % 按z方向粒子总数循环
        for j=1:size_y          % 按y方向粒子总数循环
            for k=1:size_x      % 按x方向粒子总数循环
                ii=ii+1;        % ii取值+1，从而遍历每个粒子
                particle(ii,1)=xx_min+(k-0.5)*dettx-size_boundary_left*m*dettx;     % 第ii个粒子的x方向坐标=x方向最小坐标*(k-0.5)*每个粒子长度，如果左边界定义了位移边界，第ii个粒子的x方向坐标=x方向最小坐标*(k-0.5)*每个粒子长度-粒子作用的层数*每个粒子长度
                particle(ii,2)=yy_min+(j-0.5)*dettx-size_boundary_bottom*m*dettx;	% 第ii个粒子的y方向坐标=y方向最小坐标*(j-0.5)*每个粒子长度，如果下边界定义了位移边界，第ii个粒子的y方向坐标=y方向最小坐标*(j-0.5)*每个粒子长度-粒子作用的层数*每个粒子长度
                particle(ii,3)=zz_min+(i-0.5)*dettx-size_boundary_behind*m*dettx;	% 第ii个粒子的z方向坐标=z方向最小坐标*(i-0.5)*每个粒子长度，如果后边界定义了位移边界，第ii个粒子的z方向坐标=z方向最小坐标*(i-0.5)*每个粒子长度-粒子作用的层数*每个粒子长度
            end
        end
    end
    size_particle=size(particle,1);          % 粒子总数

%%%%%%%%%%%%%%%%%%%%%
% 施加位移边界条件的粒子
% 左边界
    ii_displacement_boundary_left=0; % 左边界施加位移边界条件粒子计数赋初值
    if boundary_condition_left(1)==1 || boundary_condition_left(2)==1 || boundary_condition_left(3)==1        % 如果左边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle                                                                                 % 对所有粒子进行循环
            xx_i=particle(i,1);                                                                               % xx_i表示粒子i在x方向的坐标
            if xx_i<xx_min                                                                                    % 如果xx_i小于x方向最小坐标，说明粒子i处于左边界，则执行条件语句中的内容
                ii_displacement_boundary_left=ii_displacement_boundary_left+1;                                % 左边界施加位移边界条件粒子计数加1
                particle_boundary_left(ii_displacement_boundary_left)=i;                                      % 在集合particle_boundary_left中记录粒子i的标号
            end
        end
        size_displacement_boundary_left=ii_displacement_boundary_left;                                        % 记录左边界施加位移边界条件的粒子总数
    end
% 右边界
    ii_displacement_boundary_right=0;   % 右边界施加位移边界条件粒子计数赋初值                                
    if boundary_condition_right(1)==1 || boundary_condition_right(2)==1 || boundary_condition_right(3)==1     % 如果右边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle                                                                                 % 对所有粒子进行循环
            xx_i=particle(i,1);                                                                               % xx_i表示粒子i在x方向的坐标
            if xx_i>xx_max                                                                                    % 如果xx_i大于x方向最大坐标，说明粒子i处于右边界，则执行条件语句中的内容
                ii_displacement_boundary_right=ii_displacement_boundary_right+1;                              % 右边界施加位移边界条件粒子计数加1
                particle_boundary_right(ii_displacement_boundary_right)=i;                                    % 在集合particle_boundary_right中记录粒子i的标号
            end
        end
        size_displacement_boundary_right=ii_displacement_boundary_right;                                      % 记录右边界施加位移边界条件的粒子总数
    end
% 下边界
    ii_displacement_boundary_bottom=0;   % 下边界施加位移边界条件粒子计数赋初值
    if boundary_condition_bottom(1)==1 || boundary_condition_bottom(2)==1 || boundary_condition_bottom(3)==1  % 如果下边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle                                                                                 % 对所有粒子进行循环
            yy_i=particle(i,2);                                                                               % yy_i表示粒子i在y方向的坐标
            if yy_i<yy_min                                                                                    % 如果yy_i小于y方向最小坐标，说明粒子i处于下边界，则执行条件语句中的内容
                ii_displacement_boundary_bottom=ii_displacement_boundary_bottom+1;                            % 下边界施加位移边界条件粒子计数加1
                particle_boundary_bottom(ii_displacement_boundary_bottom)=i;                                  % 在集合particle_boundary_bottom中记录粒子i的标号
            end
        end
        size_displacement_boundary_bottom=ii_displacement_boundary_bottom;                                    % 记录下边界施加位移边界条件的粒子总数
    end
% 上边界
    ii_displacement_boundary_top=0;   % 上边界施加位移边界条件粒子计数赋初值
    if boundary_condition_top(1)==1 || boundary_condition_top(2)==1 || boundary_condition_top(3)==1           % 如果上边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle                                                                                 % 对所有粒子进行循环
            yy_i=particle(i,2);                                                                               % yy_i表示粒子i在y方向的坐标
            if yy_i>yy_max                                                                                    % 如果yy_i大于y方向最大坐标，说明粒子i处于上边界，则执行条件语句中的内容
                ii_displacement_boundary_top=ii_displacement_boundary_top+1;                                  % 上边界施加位移边界条件粒子计数加1
                particle_boundary_top(ii_displacement_boundary_top)=i;                                        % 在集合particle_boundary_top中记录粒子i的标号
            end
        end
        size_displacement_boundary_top=ii_displacement_boundary_top;                                          % 记录上边界施加均布荷载条件的粒子总数
    end
% 后边界
    ii_displacement_boundary_behind=0;   % 下边界施加位移边界条件粒子计数赋初值
    if boundary_condition_bottom(1)==1 || boundary_condition_bottom(2)==1 || boundary_condition_bottom(3)==1  % 如果后边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle                                                                                 % 对所有粒子进行循环
            zz_i=particle(i,3);                                                                               % zz_i表示粒子i在z方向的坐标
            if zz_i<zz_min                                                                                    % 如果zz_i小于z方向最小坐标，说明粒子i处于后边界，则执行条件语句中的内容
                ii_displacement_boundary_behind=ii_displacement_boundary_behind+1;                            % 后边界施加位移边界条件粒子计数加1
                particle_boundary_behind(ii_displacement_boundary_behind)=i;                                  % 在集合particle_boundary_behind中记录粒子i的标号
            end
        end
        size_displacement_boundary_behind=ii_displacement_boundary_behind;                                    % 记录后边界施加均布荷载条件的粒子总数
    end
% 前边界
    ii_displacement_boundary_front=0;   % 上边界施加位移边界条件粒子计数赋初值
    if boundary_condition_front(1)==1 || boundary_condition_front(2)==1 || boundary_condition_front(3)==1     % 如果前边界施加位移条件，则执行条件语句中的内容
        for i=1:size_particle                                                                                 % 对所有粒子进行循环
            zz_i=particle(i,3);                                                                               % zz_i表示粒子i在z方向的坐标
            if zz_i>zz_max                                                                                    % 如果zz_i大于z方向最大坐标，说明粒子i处于前边界，则执行条件语句中的内容
                ii_displacement_boundary_front=ii_displacement_boundary_front+1;                              % 前边界施加位移边界条件粒子计数加1
                particle_boundary_front(ii_displacement_boundary_front)=i;                                    % 在集合particle_boundary_front中记录粒子i的标号
            end
        end
        size_displacement_boundary_front=ii_displacement_boundary_front;                                      % 记录后前界施加均布荷载条件的粒子总数
    end
%%%%%%%%%%%%%%%%%%%%%
% 施加均布荷载条件的粒子
% 左边界
    ii_uniform_load_left=0; % 左边界施加位移边界条件粒子计数赋初值
    if boundary_condition_left(1)==2 || boundary_condition_left(2)==2 || boundary_condition_left(3)==2       % 如果左边界施加均布荷载条件，则执行条件语句中的内容
        for i=1:size_particle                                                                                % 对所有粒子进行循环
            xx_i=particle(i,1);                                                                              % xx_i表示粒子i在x方向的坐标
            if xx_i<(xx_min+dettx) && xx_i>xx_min                                                                           % 如果xx_i小于(x方向最小坐标+dettx)，说明粒子i为位于最左侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_left=ii_uniform_load_left+1;                                                 % 左边界施加均布荷载条件粒子计数加1
                particle_uniform_load_left(ii_uniform_load_left)=i;                                          % 在集合particle_uniform_load_left中记录粒子i的标号
            end
        end
        size_uniform_load_left=ii_uniform_load_left;                                                         % 记录左边界施加均布荷载条件的粒子总数
    end
% 右边界
    ii_uniform_load_right=0;   % 右边界施加位移边界条件粒子计数赋初值
    if boundary_condition_right(1)==2 || boundary_condition_right(2)==2 || boundary_condition_right(3)==2    % 如果右边界施加均布荷载条件，则执行条件语句中的内容 
        for i=1:size_particle                                                                                % 对所有粒子进行循环
            xx_i=particle(i,1);                                                                              % xx_i表示粒子i在x方向的坐标
            if xx_i>(xx_max-dettx) && xx_i<xx_max                                                                           % 如果xx_i大于(x方向最大坐标-dettx)，说明粒子i为位于最右侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_right=ii_uniform_load_right+1;                                               % 右边界施加均布荷载条件粒子计数加1
                particle_uniform_load_right(ii_uniform_load_right)=i;                                        % 在集合particle_uniform_load_right中记录粒子i的标号
            end
        end
        size_uniform_load_right=ii_uniform_load_right;                                                       % 记录右边界施加均布荷载条件的粒子总数
    end
% 下边界
    ii_uniform_load_bottom=0;   % 下边界施加位移边界条件粒子计数赋初值
    if boundary_condition_bottom(1)==2 || boundary_condition_bottom(2)==2 || boundary_condition_bottom(3)==2 % 如果下边界施加均布荷载条件，则执行条件语句中的内容
        for i=1:size_particle                                                                                % 对所有粒子进行循环
            yy_i=particle(i,2);                                                                              % yy_i表示粒子i在y方向的坐标
            if yy_i<(yy_min+dettx) && yy_i>yy_min                                                                           % 如果yy_i小于(y方向最小坐标+dettx)，说明粒子i为位于最下侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_bottom=ii_uniform_load_bottom+1;                                             % 下边界施加均布荷载条件粒子计数加1
                particle_uniform_load_bottom(ii_displacement_boundary_bottom)=i;                             % 在集合particle_uniform_load_bottom中记录粒子i的标号
            end
        end
        size_uniform_load_bottom=ii_displacement_uniform_load_bottom;                                        % 记录下边界施加均布荷载条件的粒子总数
    end
% 上边界
    ii_uniform_load_top=0;   % 上边界施加位移边界条件粒子计数赋初值
    if boundary_condition_top(1)==2 || boundary_condition_top(2)==2 || boundary_condition_top(3)==2          % 如果上边界施加均布荷载条件，则执行条件语句中的内容 
        for i=1:size_particle                                                                                % 对所有粒子进行循环
            yy_i=particle(i,2);                                                                              % yy_i表示粒子i在y方向的坐标
            if yy_i>(yy_max-dettx) && yy_i<yy_max                                                                           % 如果yy_i大于(y方向最大坐标-dettx)，说明粒子i为位于最上侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_top=ii_uniform_load_top+1;                                                   % 上边界施加均布荷载条件粒子计数加1
                particle_uniform_load_top(ii_uniform_load_top)=i;                                            % 在集合particle_uniform_load_top中记录粒子i的标号
            end
        end
        size_uniform_load_top=ii_uniform_load_top;                                                           % 记录上边界施加均布荷载条件的粒子总数
    end
% 后边界
    ii_uniform_load_behind=0;   % 下边界施加位移边界条件粒子计数赋初值
    if boundary_condition_behind(1)==2 || boundary_condition_behind(2)==2 || boundary_condition_behind(3)==2 % 如果后边界施加均布荷载条件，则执行条件语句中的内容
        for i=1:size_particle                                                                                % 对所有粒子进行循环
            zz_i=particle(i,3);                                                                              % zz_i表示粒子i在z方向的坐标
            if zz_i<(zz_min+dettx) && zz_i>zz_min                                                                           % 如果zz_i小于(z方向最小坐标+dettx)，说明粒子i为位于最后侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_behind=ii_uniform_load_behind+1;                                             % 后边界施加均布荷载条件粒子计数加1
                particle_uniform_load_behind(ii_displacement_boundary_behind)=i;                             % 在集合particle_uniform_load_behind中记录粒子i的标号
            end
        end
        size_uniform_load_behind=ii_displacement_uniform_load_behind;                                        % 记录后边界施加均布荷载条件的粒子总数
    end
% 前边界
    ii_uniform_load_front=0;   % 上边界施加位移边界条件粒子计数赋初值
    if boundary_condition_front(1)==2 || boundary_condition_front(2)==2 || boundary_condition_front(3)==2    % 如果前边界施加均布荷载条件，则执行条件语句中的内容 
        for i=1:size_particle                                                                                % 对所有粒子进行循环
            zz_i=particle(i,3);                                                                              % zz_i表示粒子i在z方向的坐标
            if zz_i>(zz_max-dettx) && zz_i<zz_max                                                                           % 如果zz_i大于(z方向最大坐标-dettx)，说明粒子i为位于最前侧一列的粒子，则执行条件语句中的内容
                ii_uniform_load_front=ii_uniform_load_front+1;                                               % 前边界施加均布荷载条件粒子计数加1
                particle_uniform_load_front(ii_uniform_load_front)=i;                                        % 在集合particle_uniform_load_front中记录粒子i的标号
            end
        end
        size_uniform_load_front=ii_uniform_load_front;                                                       % 记录前边界施加均布荷载条件的粒子总数
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

size_temperature_boundary_left=0;  % 左边界施加温度边界条件的粒子总数赋初值
particle_boundary_left_temperature=[];          % 定义左边界施加温度边界条件的粒子集合为空集
size_temperature_boundary_right=0; % 右边界施加温度边界条件的粒子总数赋初值
particle_boundary_right_temperature=[];         % 定义右边界施加温度边界条件的粒子集合为空集
size_temperature_boundary_bottom=0; % 下边界施加温度边界条件的粒子总数赋初值
particle_boundary_bottom_temperature=[];        % 定义下边界施加温度边界条件的粒子集合为空集
size_temperature_boundary_top=0;   % 上边界施加温度边界条件的粒子总数赋初值
particle_boundary_top_temperature=[];           % 定义上边界施加温度边界条件的粒子集合为空集
size_temperature_boundary_behind=0; % 后边界施加温度边界条件的粒子总数赋初值
particle_boundary_behind_temperature=[];           % 定义后边界施加温度边界条件的粒子集合为空集
size_temperature_boundary_front=0; % 前边界施加温度边界条件的粒子总数赋初值
particle_boundary_front_temperature=[];         % 定义前边界施加温度边界条件的粒子集合为空集

size_heat_flux_left=0;       % 左边界施加热流条件的粒子总数赋初值
particle_heat_flux_left=[];  % 定义左边界施加热流条件的粒子集合为空集
size_heat_flux_right=0;      % 右边界施加热流条件的粒子总数赋初值
particle_heat_flux_right=[]; % 定义右边界施加热流条件的粒子集合为空集
size_heat_flux_bottom=0;     % 下边界施加热流条件的粒子总数赋初值
particle_heat_flux_bottom=[];% 定义下边界施加热流条件的粒子集合为空集
size_heat_flux_top=0;        % 上边界施加热流条件的粒子总数赋初值
particle_heat_flux_top=[];   % 定义上边界施加热流条件的粒子集合为空集
size_heat_flux_behind=0;     % 后边界施加热流条件的粒子总数赋初值
particle_heat_flux_behind=[];% 定义后边界施加热流条件的粒子集合为空集
size_heat_flux_front=0;      % 前边界施加热流条件的粒子总数赋初值
particle_heat_flux_front=[]; % 定义前边界施加热流条件的粒子集合为空集
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%一维情况
if dimension==1
%%%%%%%%%%%%%%%%%%%%%
% 施加温度边界条件的粒子
%左边界
    ii_temperature_boundary_left=0;    % 左边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_left_temperature==1    % 如果左边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            if xx_i<xx_min              % 如果xx_i小于x方向最小坐标，说明粒子i处于左边界，则执行条件语句中的内容
                ii_temperature_boundary_left=ii_temperature_boundary_left+1;  % 左边界施加温度边界条件粒子计数加1
                particle_boundary_left_temperature(ii_temperature_boundary_left,1)=i;        % 在集合particle_boundary_left中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    if abs(xx_i+xx_j-2*xx_min)<0.1*dettx % 如果粒子i的x坐标xx_i与粒子j的x坐标xx_j的和与2倍的左边界x坐标xx_min的差的绝对值小于0.1*dettx，则粒子i与粒子j关于左边界对称，执行条件语句中的内容
                        particle_boundary_left_temperature(ii_temperature_boundary_left,2)=j;        % 在集合particle_boundary_left中记录粒子j的标号
                        break;          % 找到与粒子i关于左边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_left=ii_temperature_boundary_left;  % 记录左边界施加温度边界条件的粒子总数
    end
%右边界
    ii_temperature_boundary_right=0;    % 右边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_right_temperature==1    % 如果右边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            if xx_i>xx_max              % 如果xx_i大于x方向最大坐标，说明粒子i处于右边界，则执行条件语句中的内容
                ii_temperature_boundary_right=ii_temperature_boundary_right+1;  % 右边界施加温度边界条件粒子计数加1
                particle_boundary_right_temperature(ii_temperature_boundary_right,1)=i;        % 在集合particle_boundary_right中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    if abs(xx_i+xx_j-2*xx_max)<0.1*dettx  % 如果粒子i的x坐标xx_i与粒子j的x坐标xx_j的和与2倍的右边界x坐标xx_max的差的绝对值小于0.1*dettx，则粒子i与粒子j关于右边界对称，执行条件语句中的内容
                        particle_boundary_right_temperature(ii_temperature_boundary_right,2)=j;        % 在集合particle_boundary_right中记录粒子j的标号
                        break;          % 找到与粒子i关于右边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_right=ii_temperature_boundary_right;  % 记录右边界施加温度边界条件的粒子总数
    end
%%%%%%%%%%%%%%%%%%%%%
% 施加热流条件的粒子
%左边界
    ii_heat_flux_left=0;                                         % 左边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_left_temperature==2                                % 如果左边界施加热流条件，则执行条件语句中的内容
        for i=1:size_particle                                       % 对所有粒子进行循环
            xx_i=particle(i,1);                                     % xx_i表示粒子i在x方向的坐标
            if xx_i<(xx_min+dettx) && xx_i>xx_min                             % 如果xx_i小于(x方向最小坐标+dettx)，说明粒子i为位于最左侧一列的粒子，则执行条件语句中的内容
                ii_heat_flux_left=ii_heat_flux_left+1;        % 左边界施加热流条件粒子计数加1
                particle_heat_flux_left(ii_heat_flux_left)=i; % 在集合particle_uniform_load_left中记录粒子i的标号
            end
        end
        size_heat_flux_left=ii_heat_flux_left;                % 记录左边界施加热流条件的粒子总数
    end
%右边界
    ii_heat_flux_right=0;                                        % 右边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_right_temperature==2                               % 如果右边界施加热流条件，则执行条件语句中的内容 
        for i=1:size_particle                                       % 对所有粒子进行循环
            xx_i=particle(i,1);                                     % xx_i表示粒子i在x方向的坐标
            if xx_i>(xx_max-dettx) && xx_i<xx_max                                   % 如果xx_i大于(x方向最大坐标-dettx)，说明粒子i为位于最右侧一列的粒子，则执行条件语句中的内容
                ii_heat_flux_right=ii_heat_flux_right+1;        % 右边界施加热流条件粒子计数加1
                particle_heat_flux_right(ii_heat_flux_right)=i; % 在集合particle_uniform_load_right中记录粒子i的标号
        	end
        end
        size_heat_flux_right=ii_heat_flux_right;                % 记录右边界施加热流条件的粒子总数
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%二维情况
if dimension==2
%%%%%%%%%%%%%%%%%%%%%
% 施加温度边界条件的粒子
%左边界
    ii_temperature_boundary_left=0;    % 左边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_left_temperature==1    % 如果左边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            yy_i=particle(i,2);         % yy_i表示粒子i在y方向的坐标
            if xx_i<xx_min              % 如果xx_i小于x方向最小坐标，说明粒子i处于左边界，则执行条件语句中的内容
                ii_temperature_boundary_left=ii_temperature_boundary_left+1;  % 左边界施加温度边界条件粒子计数加1
                particle_boundary_left_temperature(ii_temperature_boundary_left,1)=i;        % 在集合particle_boundary_left中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    yy_j=particle(j,2);         % yy_j表示粒子j在y方向的坐标
                    if abs(xx_i+xx_j-2*xx_min)<0.1*dettx && abs(yy_i-yy_j)<0.1*dettx % 如果粒子i的x坐标xx_i与粒子j的x坐标xx_j的和与2倍的左边界x坐标xx_min的差的绝对值小于0.1*dettx，且粒子i的y坐标yy_i与粒子j的y坐标yy_j的差的绝对值小于0.1*dettx，则粒子i与粒子j关于左边界对称，执行条件语句中的内容
                        particle_boundary_left_temperature(ii_temperature_boundary_left,2)=j;        % 在集合particle_boundary_left中记录粒子j的标号
                        break;          % 找到与粒子i关于左边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_left=ii_temperature_boundary_left;  % 记录左边界施加温度边界条件的粒子总数
    end
%右边界
    ii_temperature_boundary_right=0;    % 右边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_right_temperature==1    % 如果右边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            yy_i=particle(i,2);         % yy_i表示粒子i在y方向的坐标
            if xx_i>xx_max              % 如果xx_i大于x方向最大坐标，说明粒子i处于右边界，则执行条件语句中的内容
                ii_temperature_boundary_right=ii_temperature_boundary_right+1;  % 右边界施加温度边界条件粒子计数加1
                particle_boundary_right_temperature(ii_temperature_boundary_right,1)=i;        % 在集合particle_boundary_right中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    yy_j=particle(j,2);         % yy_j表示粒子j在y方向的坐标
                    if abs(xx_i+xx_j-2*xx_max)<0.1*dettx && abs(yy_i-yy_j)<0.1*dettx % 如果粒子i的x坐标xx_i与粒子j的x坐标xx_j的和与2倍的右边界x坐标xx_max的差的绝对值小于0.1*dettx，且粒子i的y坐标yy_i与粒子j的y坐标yy_j的差的绝对值小于0.1*dettx，则粒子i与粒子j关于右边界对称，执行条件语句中的内容
                        particle_boundary_right_temperature(ii_temperature_boundary_right,2)=j;        % 在集合particle_boundary_right中记录粒子j的标号
                        break;          % 找到与粒子i关于右边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_right=ii_temperature_boundary_right;  % 记录右边界施加温度边界条件的粒子总数
    end
%下边界
    ii_temperature_boundary_bottom=0;   % 下边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_bottom_temperature==1     % 如果下边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            yy_i=particle(i,2);         % yy_i表示粒子i在y方向的坐标
            if yy_i<yy_min              % 如果yy_i小于y方向最小坐标，说明粒子i处于下边界，则执行条件语句中的内容
                ii_temperature_boundary_bottom=ii_temperature_boundary_bottom+1;  % 下边界施加温度边界条件粒子计数加1
                particle_boundary_bottom_temperature(ii_temperature_boundary_bottom,1)=i;        % 在集合particle_boundary_bottom中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    yy_j=particle(j,2);         % yy_j表示粒子j在y方向的坐标
                    if abs(yy_i+yy_j-2*yy_min)<0.1*dettx && abs(xx_i-xx_j)<0.1*dettx % 如果粒子i的y坐标yy_i与粒子j的y坐标yy_j的和与2倍的下边界y坐标yy_min的差的绝对值小于0.1*dettx，且粒子i的x坐标xx_i与粒子j的x坐标xx_j的差的绝对值小于0.1*dettx，则粒子i与粒子j关于下边界对称，执行条件语句中的内容
                        particle_boundary_bottom_temperature(ii_temperature_boundary_bottom,2)=j;        % 在集合particle_boundary_bottom中记录粒子j的标号
                        break;          % 找到与粒子i关于下边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_bottom=ii_temperature_boundary_bottom;  % 记录下边界施加温度边界条件的粒子总数
    end
%上边界
    ii_temperature_boundary_top=0;   % 上边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_top_temperature==1     % 如果上边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            yy_i=particle(i,2);         % yy_i表示粒子i在y方向的坐标
            if yy_i>yy_max              % 如果yy_i大于y方向最大坐标，说明粒子i处于上边界，则执行条件语句中的内容
                ii_temperature_boundary_top=ii_temperature_boundary_top+1;  % 上边界施加温度边界条件粒子计数加1
                particle_boundary_top_temperature(ii_temperature_boundary_top,1)=i;       % 在集合particle_boundary_top中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    yy_j=particle(j,2);         % yy_j表示粒子j在y方向的坐标
                    if abs(yy_i+yy_j-2*yy_max)<0.1*dettx && abs(xx_i-xx_j)<0.1*dettx % 如果粒子i的y坐标yy_i与粒子j的y坐标yy_j的和与2倍的上边界y坐标yy_max的差的绝对值小于0.1*dettx，且粒子i的x坐标xx_i与粒子j的x坐标xx_j的差的绝对值小于0.1*dettx，则粒子i与粒子j关于上边界对称，执行条件语句中的内容
                        particle_boundary_top_temperature(ii_temperature_boundary_top,2)=j;        % 在集合particle_boundary_top中记录粒子j的标号
                        break;          % 找到与粒子i关于上边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_top=ii_temperature_boundary_top;  % 记录上边界施加温度边界条件的粒子总数
    end
%%%%%%%%%%%%%%%%%%%%%
% 施加热流条件的粒子
%左边界
    ii_heat_flux_left=0;                                         % 左边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_left_temperature==2                                % 如果左边界施加热流条件，则执行条件语句中的内容
        for i=1:size_particle                                       % 对所有粒子进行循环
            xx_i=particle(i,1);                                     % xx_i表示粒子i在x方向的坐标
            if xx_i<(xx_min+dettx) && xx_i>xx_min                                   % 如果xx_i小于(x方向最小坐标+dettx)，说明粒子i为位于最左侧一列的粒子，则执行条件语句中的内容
                ii_heat_flux_left=ii_heat_flux_left+1;        % 左边界施加热流条件粒子计数加1
                particle_heat_flux_left(ii_heat_flux_left)=i; % 在集合particle_uniform_load_left中记录粒子i的标号
            end
        end
        size_heat_flux_left=ii_heat_flux_left;                % 记录左边界施加热流条件的粒子总数
    end
%右边界
    ii_heat_flux_right=0;                                        % 右边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_right_temperature==2                               % 如果右边界施加热流条件，则执行条件语句中的内容 
        for i=1:size_particle                                       % 对所有粒子进行循环
            xx_i=particle(i,1);                                     % xx_i表示粒子i在x方向的坐标
            if xx_i>(xx_max-dettx) && xx_i<xx_max                                  % 如果xx_i大于(x方向最大坐标-dettx)，说明粒子i为位于最右侧一列的粒子，则执行条件语句中的内容
                ii_heat_flux_right=ii_heat_flux_right+1;        % 右边界施加热流条件粒子计数加1
                particle_heat_flux_right(ii_heat_flux_right)=i; % 在集合particle_uniform_load_right中记录粒子i的标号
        	end
        end
        size_heat_flux_right=ii_heat_flux_right;                % 记录右边界施加热流条件的粒子总数
    end

%下边界
    ii_heat_flux_bottom=0;                                        % 下边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_bottom_temperature==2                               % 如果下边界施加热流条件，则执行条件语句中的内容 
        for i=1:size_particle                                       % 对所有粒子进行循环
            yy_i=particle(i,2);                                     % yy_i表示粒子i在y方向的坐标
            if yy_i<(yy_min+dettx) && yy_i>yy_min                                  % 如果yy_i小于(y方向最小坐标+dettx)，说明粒子i为位于最下侧一行的粒子，则执行条件语句中的内容
                ii_heat_flux_bottom=ii_heat_flux_bottom+1;        % 下边界施加热流条件粒子计数加1
                particle_heat_flux_bottom(ii_heat_flux_bottom)=i; % 在集合particle_uniform_load_bottom中记录粒子i的标号
            end
        end
        size_heat_flux_bottom=ii_heat_flux_bottom;                % 记录下边界施加热流条件的粒子总数
    end
%上边界
    ii_heat_flux_top=0;                                        % 上边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_top_temperature==2                               % 如果上边界施加热流条件，则执行条件语句中的内容
        for i=1:size_particle                                       % 对所有粒子进行循环
            yy_i=particle(i,2);                                     % yy_i表示粒子i在y方向的坐标
            if yy_i>(yy_max-dettx) && yy_i<yy_max                                 % 如果yy_i大于(y方向最大坐标-dettx)，说明粒子i为位于最上侧一行的粒子，则执行条件语句中的内容
                ii_heat_flux_top=ii_heat_flux_top+1;        % 上边界施加热流条件粒子计数加1
                particle_heat_flux_top(ii_heat_flux_top)=i; % 在集合particle_uniform_load_top中记录粒子i的标号
            end
        end
        size_heat_flux_top=ii_heat_flux_top;                % 记录上边界施加热流条件的粒子总数
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%三维情况
if dimension==3
%%%%%%%%%%%%%%%%%%%%%
% 施加温度边界条件的粒子
%左边界
    ii_temperature_boundary_left=0;    % 左边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_left_temperature==1    % 如果左边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            yy_i=particle(i,2);         % yy_i表示粒子i在y方向的坐标
            zz_i=particle(i,3);         % zz_i表示粒子i在z方向的坐标
            if xx_i<xx_min              % 如果xx_i小于x方向最小坐标，说明粒子i处于左边界，则执行条件语句中的内容
                ii_temperature_boundary_left=ii_temperature_boundary_left+1;  % 左边界施加温度边界条件粒子计数加1
                particle_boundary_left_temperature(ii_temperature_boundary_left,1)=i;        % 在集合particle_boundary_left中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    yy_j=particle(j,2);         % yy_j表示粒子j在y方向的坐标
                    zz_j=particle(j,3);         % zz_j表示粒子j在z方向的坐标
                    if abs(xx_i+xx_j-2*xx_min)<0.1*dettx && abs(yy_i-yy_j)<0.1*dettx && abs(zz_i-zz_j)<0.1*dettx % 如果粒子i的x坐标xx_i与粒子j的x坐标xx_j的和与2倍的左边界x坐标xx_min的差的绝对值小于0.1*dettx，且粒子i的y坐标yy_i与粒子j的y坐标yy_j的差的绝对值小于0.1*dettx，且粒子i的z坐标zz_i与粒子j的z坐标zz_j的差的绝对值小于0.1*dettx，则粒子i与粒子j关于左边界对称，执行条件语句中的内容
                        particle_boundary_left_temperature(ii_temperature_boundary_left,2)=j;        % 在集合particle_boundary_left中记录粒子j的标号
                        break;          % 找到与粒子i关于左边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_left=ii_temperature_boundary_left;  % 记录左边界施加温度边界条件的粒子总数
    end
%右边界
    ii_temperature_boundary_right=0;    % 右边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_right_temperature==1    % 如果右边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            yy_i=particle(i,2);         % yy_i表示粒子i在y方向的坐标
            zz_i=particle(i,3);         % zz_i表示粒子i在z方向的坐标
            if xx_i>xx_max              % 如果xx_i大于x方向最大坐标，说明粒子i处于右边界，则执行条件语句中的内容
                ii_temperature_boundary_right=ii_temperature_boundary_right+1;  % 右边界施加温度边界条件粒子计数加1
                particle_boundary_right_temperature(ii_temperature_boundary_right,1)=i;        % 在集合particle_boundary_right中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    yy_j=particle(j,2);         % yy_j表示粒子j在y方向的坐标
                    zz_j=particle(j,3);         % zz_j表示粒子j在z方向的坐标
                    if abs(xx_i+xx_j-2*xx_max)<0.1*dettx && abs(yy_i-yy_j)<0.1*dettx && abs(zz_i-zz_j)<0.1*dettx % 如果粒子i的x坐标xx_i与粒子j的x坐标xx_j的和与2倍的右边界x坐标xx_max的差的绝对值小于0.1*dettx，且粒子i的y坐标yy_i与粒子j的y坐标yy_j的差的绝对值小于0.1*dettx，且粒子i的z坐标zz_i与粒子j的z坐标zz_j的差的绝对值小于0.1*dettx，则粒子i与粒子j关于右边界对称，执行条件语句中的内容
                        particle_boundary_right_temperature(ii_temperature_boundary_right,2)=j;        % 在集合particle_boundary_right中记录粒子j的标号
                        break;          % 找到与粒子i关于右边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_right=ii_temperature_boundary_right;  % 记录右边界施加温度边界条件的粒子总数
    end
%下边界
    ii_temperature_boundary_bottom=0;   % 下边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_bottom_temperature==1     % 如果下边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            yy_i=particle(i,2);         % yy_i表示粒子i在y方向的坐标
            zz_i=particle(i,3);         % zz_i表示粒子i在z方向的坐标
            if yy_i<yy_min              % 如果yy_i小于y方向最小坐标，说明粒子i处于下边界，则执行条件语句中的内容
                ii_temperature_boundary_bottom=ii_temperature_boundary_bottom+1;  % 下边界施加温度边界条件粒子计数加1
                particle_boundary_bottom_temperature(ii_temperature_boundary_bottom,1)=i;        % 在集合particle_boundary_bottom中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    yy_j=particle(j,2);         % yy_j表示粒子j在y方向的坐标
                    zz_j=particle(j,3);         % zz_j表示粒子j在z方向的坐标
                    if abs(yy_i+yy_j-2*yy_min)<0.1*dettx && abs(xx_i-xx_j)<0.1*dettx && abs(zz_i-zz_j)<0.1*dettx % 如果粒子i的y坐标yy_i与粒子j的y坐标yy_j的和与2倍的下边界y坐标yy_min的差的绝对值小于0.1*dettx，且粒子i的x坐标xx_i与粒子j的x坐标xx_j的差的绝对值小于0.1*dettx，且粒子i的z坐标zz_i与粒子j的z坐标zz_j的差的绝对值小于0.1*dettx,则粒子i与粒子j关于下边界对称，执行条件语句中的内容
                        particle_boundary_bottom_temperature(ii_temperature_boundary_bottom,2)=j;        % 在集合particle_boundary_bottom中记录粒子j的标号
                        break;          % 找到与粒子i关于下边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_bottom=ii_temperature_boundary_bottom;  % 记录下边界施加温度边界条件的粒子总数
    end
%上边界
    ii_temperature_boundary_top=0;   % 上边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_top_temperature==1     % 如果上边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            yy_i=particle(i,2);         % yy_i表示粒子i在y方向的坐标
            zz_i=particle(i,3);         % zz_i表示粒子i在z方向的坐标
            if yy_i>yy_max              % 如果yy_i大于y方向最大坐标，说明粒子i处于上边界，则执行条件语句中的内容
                ii_temperature_boundary_top=ii_temperature_boundary_top+1;  % 上边界施加温度边界条件粒子计数加1
                particle_boundary_top_temperature(ii_temperature_boundary_top,1)=i;       % 在集合particle_boundary_top中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    yy_j=particle(j,2);         % yy_j表示粒子j在y方向的坐标
                    zz_j=particle(j,3);         % zz_j表示粒子j在z方向的坐标
                    if abs(yy_i+yy_j-2*yy_max)<0.1*dettx && abs(xx_i-xx_j)<0.1*dettx && abs(zz_i-zz_j)<0.1*dettx % 如果粒子i的y坐标yy_i与粒子j的y坐标yy_j的和与2倍的上边界y坐标yy_max的差的绝对值小于0.1*dettx，且粒子i的x坐标xx_i与粒子j的x坐标xx_j的差的绝对值小于0.1*dettx，且粒子i的z坐标zz_i与粒子j的z坐标zz_j的差的绝对值小于0.1*dettx，则粒子i与粒子j关于上边界对称，执行条件语句中的内容
                        particle_boundary_top_temperature(ii_temperature_boundary_top,2)=j;        % 在集合particle_boundary_top中记录粒子j的标号
                        break;          % 找到与粒子i关于上边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_top=ii_temperature_boundary_top;  % 记录上边界施加温度边界条件的粒子总数
    end
% 后边界
    ii_temperature_boundary_behind=0;   % 后边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_behind_temperature==1     % 如果后边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            yy_i=particle(i,2);         % yy_i表示粒子i在y方向的坐标
            zz_i=particle(i,3);         % zz_i表示粒子i在z方向的坐标
            if zz_i<zz_min              % 如果zz_i小于z方向最小坐标，说明粒子i处于后边界，则执行条件语句中的内容
                ii_temperature_boundary_behind=ii_temperature_boundary_behind+1;  % 后边界施加温度边界条件粒子计数加1
                particle_boundary_behind_temperature(ii_temperature_boundary_behind,1)=i;       % 在集合particle_boundary_behind中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    yy_j=particle(j,2);         % yy_j表示粒子j在y方向的坐标
                    zz_j=particle(j,3);         % zz_j表示粒子j在z方向的坐标
                    if abs(zz_i+zz_j-2*zz_min)<0.1*dettx && abs(xx_i-xx_j)<0.1*dettx && abs(yy_i-yy_j)<0.1*dettx % 如果粒子i的z坐标zz_i与粒子j的z坐标zz_j的和与2倍的后边界z坐标zz_min的差的绝对值小于0.1*dettx，且粒子i的x坐标xx_i与粒子j的x坐标xx_j的差的绝对值小于0.1*dettx，且粒子i的y坐标yy_i与粒子j的y坐标yy_j的差的绝对值小于0.1*dettx,则粒子i与粒子j关于后边界对称，执行条件语句中的内容
                        particle_boundary_behind_temperature(ii_temperature_boundary_behind,2)=j;        % 在集合particle_boundary_behind中记录粒子j的标号
                        break;          % 找到与粒子i关于后边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_behind=ii_temperature_boundary_behind;  % 记录后边界施加温度边界条件的粒子总数
    end
% 前边界
    ii_temperature_boundary_front=0;   % 前边界施加温度边界条件的粒子计数赋初值为0
    if boundary_condition_front_temperature==1     % 如果前边界施加温度边界条件，则执行条件语句中的内容
        for i=1:size_particle           % 对所有粒子进行循环
            xx_i=particle(i,1);         % xx_i表示粒子i在x方向的坐标
            yy_i=particle(i,2);         % yy_i表示粒子i在y方向的坐标
            zz_i=particle(i,3);         % zz_i表示粒子i在z方向的坐标
            if zz_i>zz_max              % 如果zz_i大于z方向最大坐标，说明粒子i处于前边界，则执行条件语句中的内容
                ii_temperature_boundary_front=ii_temperature_boundary_front+1;  % 前边界施加温度边界条件粒子计数加1
                particle_boundary_front_temperature(ii_temperature_boundary_front,1)=i;       % 在集合particle_boundary_front中记录粒子i的标号
                for j=1:size_particle           % 对所有粒子进行循环
                    xx_j=particle(j,1);         % xx_j表示粒子j在x方向的坐标
                    yy_j=particle(j,2);         % yy_j表示粒子j在y方向的坐标
                    zz_j=particle(j,3);         % zz_j表示粒子j在z方向的坐标
                    if abs(zz_i+zz_j-2*zz_max)<0.1*dettx && abs(xx_i-xx_j)<0.1*dettx && abs(yy_i-yy_j)<0.1*dettx % 如果粒子i的z坐标zz_i与粒子j的z坐标zz_j的和与2倍的前边界z坐标zz_max的差的绝对值小于0.1*dettx，且粒子i的x坐标xx_i与粒子j的x坐标xx_j的差的绝对值小于0.1*dettx，且粒子i的y坐标yy_i与粒子j的y坐标yy_j的差的绝对值小于0.1*dettx,则粒子i与粒子j关于前边界对称，执行条件语句中的内容
                        particle_boundary_front_temperature(ii_temperature_boundary_front,2)=j;        % 在集合particle_boundary_front中记录粒子j的标号
                        break;          % 找到与粒子i关于前边界对称的粒子j后，跳出该循环
                    end
                end
            end
        end
        size_temperature_boundary_front=ii_temperature_boundary_front;  % 记录前边界施加温度边界条件的粒子总数
    end
%%%%%%%%%%%%%%%%%%%%%
% 施加均布荷载条件的粒子
%左边界
    ii_heat_flux_left=0;                                         % 左边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_left_temperature==2                                % 如果左边界施加热流条件，则执行条件语句中的内容
        for i=1:size_particle                                       % 对所有粒子进行循环
            xx_i=particle(i,1);                                     % xx_i表示粒子i在x方向的坐标
            if xx_i<(xx_min+dettx) && xx_i>xx_min                                 % 如果xx_i小于(x方向最小坐标+dettx)，说明粒子i为位于最左侧一面的粒子，则执行条件语句中的内容
                ii_heat_flux_left=ii_heat_flux_left+1;        % 左边界施加热流条件粒子计数加1
                particle_heat_flux_left(ii_heat_flux_left)=i; % 在集合particle_uniform_load_left中记录粒子i的标号
            end
        end
        size_heat_flux_left=ii_heat_flux_left;                % 记录左边界施加热流条件的粒子总数
    end
%右边界
    ii_heat_flux_right=0;                                        % 右边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_right_temperature==2                               % 如果右边界施加热流条件，则执行条件语句中的内容 
        for i=1:size_particle                                       % 对所有粒子进行循环
            xx_i=particle(i,1);                                     % xx_i表示粒子i在x方向的坐标
            if xx_i>(xx_max-dettx) && xx_i<xx_max                                  % 如果xx_i大于(x方向最大坐标-dettx)，说明粒子i为位于最右侧一面的粒子，则执行条件语句中的内容
                ii_heat_flux_right=ii_heat_flux_right+1;        % 右边界施加热流条件粒子计数加1
                particle_heat_flux_right(ii_heat_flux_right)=i; % 在集合particle_uniform_load_right中记录粒子i的标号
        	end
        end
        size_heat_flux_right=ii_heat_flux_right;                % 记录右边界施加热流条件的粒子总数
    end

%下边界
    ii_heat_flux_bottom=0;                                        % 下边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_bottom_temperature==2                               % 如果下边界施加热流条件，则执行条件语句中的内容 
        for i=1:size_particle                                       % 对所有粒子进行循环
            yy_i=particle(i,2);                                     % yy_i表示粒子i在y方向的坐标
            if yy_i<(yy_min+dettx) && yy_i>yy_min                                  % 如果yy_i小于(y方向最小坐标+dettx)，说明粒子i为位于最下侧一面的粒子，则执行条件语句中的内容
                ii_heat_flux_bottom=ii_heat_flux_bottom+1;        % 下边界施加热流条件粒子计数加1
                particle_heat_flux_bottom(ii_heat_flux_bottom)=i; % 在集合particle_uniform_load_bottom中记录粒子i的标号
            end
        end
        size_heat_flux_bottom=ii_heat_flux_bottom;                % 记录下边界施加热流条件的粒子总数
    end
%上边界
    ii_heat_flux_top=0;                                        % 上边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_top_temperature==2                               % 如果上边界施加热流条件，则执行条件语句中的内容
        for i=1:size_particle                                       % 对所有粒子进行循环
            yy_i=particle(i,2);                                     % yy_i表示粒子i在y方向的坐标
            if yy_i>(yy_max-dettx) && yy_i<yy_max                                 % 如果yy_i大于(y方向最大坐标-dettx)，说明粒子i为位于最上侧一面的粒子，则执行条件语句中的内容
                ii_heat_flux_top=ii_heat_flux_top+1;        % 上边界施加热流条件粒子计数加1
                particle_heat_flux_top(ii_heat_flux_top)=i; % 在集合particle_uniform_load_top中记录粒子i的标号
            end
        end
        size_heat_flux_top=ii_heat_flux_top;                % 记录上边界施加热流条件的粒子总数
    end
% 后边界
    ii_heat_flux_behind=0;                                        % 后边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_behind_temperature==2                               % 如果后边界施加热流条件，则执行条件语句中的内容
        for i=1:size_particle                                       % 对所有粒子进行循环
            zz_i=particle(i,3);                                     % zz_i表示粒子i在z方向的坐标
            if zz_i<(zz_min+dettx) && zz_i>zz_min                                  % 如果zz_i小于(z方向最小坐标+dettx)，说明粒子i为位于最后侧一面的粒子，则执行条件语句中的内容
                ii_heat_flux_behind=ii_heat_flux_behind+1;          % 后边界施加热流条件粒子计数加1
                particle_heat_flux_behind(ii_heat_flux_behind)=i;   % 在集合particle_uniform_load_behind中记录粒子i的标号
            end
        end
        size_heat_flux_behind=ii_heat_flux_behind;                % 记录后边界施加热流条件的粒子总数
    end
% 前边界
    ii_heat_flux_front=0;                                        % 前边界施加热流条件的粒子计数赋初值为0
    if boundary_condition_front_temperature==2                               % 如果前边界施加热流条件，则执行条件语句中的内容
        for i=1:size_particle                                       % 对所有粒子进行循环
            zz_i=particle(i,3);                                     % zz_i表示粒子i在z方向的坐标
            if zz_i>(zz_max-dettx) && zz_i<zz_max                                  % 如果zz_i大于(z方向最大坐标-dettx)，说明粒子i为位于最前侧一面的粒子，则执行条件语句中的内容
                ii_heat_flux_front=ii_heat_flux_front+1;          % 前边界施加热流条件粒子计数加1
                particle_heat_flux_front(ii_heat_flux_front)=i;   % 在集合particle_uniform_load_front中记录粒子i的标号
            end
        end
        size_heat_flux_front=ii_heat_flux_front;                % 记录前边界施加热流条件的粒子总数
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global_velocity_vector=zeros(size_particle*dimension,1);    % 总体速度矢量赋初值为0
global_displacement_vector=zeros(size_particle*dimension,1);    % 总体位移矢量赋初值为0
global_temperature_vector=zeros(size_particle,1);    % 总体温度矢量赋初值为0
for i=1:size_particle                                       % 对所有粒子进行循环
    global_temperature_vector(i,1)=temperature_initial; % 为粒子i赋予初始温度
end