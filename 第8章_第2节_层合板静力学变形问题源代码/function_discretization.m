function [particle,size_particle,size_x,size_y,size_z,size_displacement_boundary_left,particle_boundary_left,size_displacement_boundary_right,particle_boundary_right, ...
          size_displacement_boundary_bottom,particle_boundary_bottom,size_displacement_boundary_top,particle_boundary_top,size_displacement_boundary_behind, ...
          particle_boundary_behind,size_displacement_boundary_front,particle_boundary_front,size_uniform_load_left,particle_uniform_load_left,size_uniform_load_right, ...
          particle_uniform_load_right,size_uniform_load_bottom,particle_uniform_load_bottom,size_uniform_load_top,particle_uniform_load_top,size_uniform_load_behind, ...
          particle_uniform_load_behind,size_uniform_load_front,particle_uniform_load_front,particle_ply]=function_discretization(boundary_condition_left, ...
          boundary_condition_right,boundary_condition_bottom,boundary_condition_top,xx_min,yy_min,zz_min,xx_max,yy_max,length_x,length_y,dettx,m,thick,size_ply)

size_displacement_boundary_left=0;  % 左边界施加位移边界条件的粒子总数赋初值
particle_boundary_left=[];          % 定义左边界施加位移边界条件的粒子集合为空集0
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
size_boundary_x=0;                      % 用来统计x方向施加位移边界条件的个数，该变量初始值取0
size_boundary_left=0;                   % 用于表示左边界是否施加位移条件，=0表示左边界未施加位移条件，=1表示左边界施加位移条件，该变量初始值取0
if boundary_condition_left(1)==1 || boundary_condition_left(2)==1 || boundary_condition_left(3)==1                % 如果左边界施加位移条件，则执行条件语句中的内容
	size_boundary_x=size_boundary_x+1;      % x方向施加位移边界条件的个数+1
	size_boundary_left=1;                   % =1 代表左边界施加位移条件
end
if boundary_condition_right(1)==1 || boundary_condition_right(2)==1 || boundary_condition_right(3)==1               % 如果右边界施加位移条件，则执行条件语句中的内容
	size_boundary_x=size_boundary_x+1;      % x方向施加位移边界条件的个数+1
end
size_boundary_y=0;                      % 用来统计y方向施加位移边界条件的个数，该变量初始值取0
size_boundary_bottom=0;                 % 用于表示下边界是否施加位移条件，=0表示下边界未施加位移条件，=1表示下边界施加位移条件，该变量初始值取0
if boundary_condition_bottom(1)==1 || boundary_condition_bottom(2)==1 || boundary_condition_bottom(3)==1         % 如果下边界施加位移条件，则执行条件语句中的内容
    size_boundary_y=size_boundary_y+1;  % y方向施加位移边界条件的个数+1
	size_boundary_bottom=1;              % =1 代表下边界施加位移条件
end
if boundary_condition_top(1)==1 || boundary_condition_top(2)==1 || boundary_condition_top(3)==1            % 如果上边界施加位移条件，则执行条件语句中的内容
	size_boundary_y=size_boundary_y+1;  % y方向施加位移边界条件的个数+1
end

size_x=length_x/dettx+size_boundary_x*m;     % x方向粒子总数=模型x方向长度/每个粒子长度+x方向施加位移边界条件的个数*粒子作用的层数
size_y=length_y/dettx+size_boundary_y*m;     % y方向粒子总数=模型y方向长度/每个粒子长度+y方向施加位移边界条件的个数*粒子作用的层数
size_z=size_ply;        % z方向粒子总数=层合板层数                                   % 
ii=0;       % 变量初始值取0
particle=zeros(size_x*size_y*size_z,3);   % 所有粒子坐标赋初值为0
particle_ply=zeros(size_x*size_y*size_z,1); % 所有粒子所在层数赋初值为0
for i=1:size_z              % 按z方向粒子总数循环
	for j=1:size_y          % 按y方向粒子总数循环
        for k=1:size_x      % 按x方向粒子总数循环
            ii=ii+1;        % ii取值+1，从而遍历每个粒子
            particle(ii,1)=xx_min+(k-0.5)*dettx-size_boundary_left*m*dettx;     % 第ii个粒子的x方向坐标=x方向最小坐标*(k-0.5)*每个粒子长度，如果左边界定义了位移边界，第ii个粒子的x方向坐标=x方向最小坐标*(k-0.5)*每个粒子长度-粒子作用的层数*每个粒子长度
            particle(ii,2)=yy_min+(j-0.5)*dettx-size_boundary_bottom*m*dettx;	% 第ii个粒子的y方向坐标=y方向最小坐标*(j-0.5)*每个粒子长度，如果下边界定义了位移边界，第ii个粒子的y方向坐标=y方向最小坐标*(j-0.5)*每个粒子长度-粒子作用的层数*每个粒子长度
            particle(ii,3)=zz_min+(i-0.5)*thick;	% 第ii个粒子的z方向坐标=z方向最小坐标*(i-0.5)*单层厚度
            particle_ply(ii,1)=i;	% 第ii个粒子属于第i层
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
%%%%%%%%%%%%%%%%%%%%%
% 施加均布荷载条件的粒子
% 左边界
ii_uniform_load_left=0; % 左边界施加位移边界条件粒子计数赋初值
if boundary_condition_left(1)==2 || boundary_condition_left(2)==2 || boundary_condition_left(3)==2       % 如果左边界施加均布荷载条件，则执行条件语句中的内容
	for i=1:size_particle                                                                                % 对所有粒子进行循环
        xx_i=particle(i,1);                                                                              % xx_i表示粒子i在x方向的坐标
        if xx_i<(xx_min+dettx)                                                                           % 如果xx_i小于(x方向最小坐标+dettx)，说明粒子i为位于最左侧一列的粒子，则执行条件语句中的内容
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
        if xx_i>(xx_max-dettx)                                                                           % 如果xx_i大于(x方向最大坐标-dettx)，说明粒子i为位于最右侧一列的粒子，则执行条件语句中的内容
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
        if yy_i<(yy_min+dettx)                                                                           % 如果yy_i小于(y方向最小坐标+dettx)，说明粒子i为位于最下侧一列的粒子，则执行条件语句中的内容
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
        if yy_i>(yy_max-dettx)                                                                           % 如果yy_i大于(y方向最大坐标-dettx)，说明粒子i为位于最上侧一列的粒子，则执行条件语句中的内容
            ii_uniform_load_top=ii_uniform_load_top+1;                                                   % 上边界施加均布荷载条件粒子计数加1
            particle_uniform_load_top(ii_uniform_load_top)=i;                                            % 在集合particle_uniform_load_top中记录粒子i的标号
        end
    end
	size_uniform_load_top=ii_uniform_load_top;                                                           % 记录上边界施加均布荷载条件的粒子总数
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%