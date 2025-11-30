function[global_displacement_vector]=function_displacement_solution(dimension,size_particle,global_stiffness_matrix,size_displacement_boundary_left,particle_boundary_left, ...
         size_displacement_boundary_right,particle_boundary_right,size_displacement_boundary_bottom,particle_boundary_bottom,size_displacement_boundary_top, ...
         particle_boundary_top,size_displacement_boundary_behind,particle_boundary_behind,size_displacement_boundary_front,particle_boundary_front, ...
         size_uniform_load_left,particle_uniform_load_left,size_uniform_load_right,particle_uniform_load_right,size_uniform_load_bottom,particle_uniform_load_bottom, ...
         size_uniform_load_top,particle_uniform_load_top,size_uniform_load_behind,particle_uniform_load_behind,size_uniform_load_front,particle_uniform_load_front, ...
         displacement_left,displacement_right,displacement_bottom,displacement_top,displacement_behind,displacement_front,uniform_load_left,uniform_load_right, ...
         uniform_load_bottom,uniform_load_top,uniform_load_behind,uniform_load_front,body_force,boundary_condition_left,boundary_condition_right,boundary_condition_bottom, ...
         boundary_condition_top,boundary_condition_behind,boundary_condition_front,volume_particle)

global_load_vector=zeros(size_particle*dimension,1); % 总体荷载向量赋初值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%一维情况
if dimension==1
%%%%%%%%%%%%%%%%%%%%%%%
%施加位移边界条件
    % 左边界
    if boundary_condition_left(1)==1                                % 如果左边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left                     % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);                           % ii表示左边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii,ii)=global_stiffness_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii,1)=global_stiffness_matrix(ii,ii)*displacement_left(1);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    % 右边界
    if boundary_condition_right(1)==1                           % 如果右边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right                % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);                      % ii表示右边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii,ii)=global_stiffness_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii,1)=global_stiffness_matrix(ii,ii)*displacement_right(1); % 采用乘大数法对总体荷载向量进行操作
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加均布荷载条件
    % 左边界
    if boundary_condition_left(1)==2    % 如果左边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left  % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);   % ii表示左边界施加荷载条件的粒子编号
            global_load_vector(ii,1)=global_load_vector(ii,1)+uniform_load_left(1)*volume_particle; % 在总体荷载向量对应位置增加荷载
        end
    end
    % 右边界
    if boundary_condition_right(1)==2    % 如果右边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right  % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);   % ii表示右边界施加荷载条件的粒子编号
            global_load_vector(ii,1)=global_load_vector(ii,1)+uniform_load_right(1)*volume_particle; % 在总体荷载向量对应位置增加荷载
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加体力
    if body_force(1)~=0    % 如果x方向体力密度不为0，说明施加了x方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle % 对所有粒子循环
            global_load_vector(ii,1)=global_load_vector(ii,1)+body_force(1)*volume_particle; % 在总体荷载向量对应位置增加荷载
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%二维情况
if dimension==2
%%%%%%%%%%%%%%%%%%%%%%%
%施加位移边界条件
    % 左边界
    if boundary_condition_left(1)==1    % 如果左边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left    % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*2-1,ii*2-1)=global_stiffness_matrix(ii*2-1,ii*2-1)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*2-1,1)=global_stiffness_matrix(ii*2-1,ii*2-1)*displacement_left(1);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_left(2)==1    % 如果左边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left    % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*2,ii*2)=global_stiffness_matrix(ii*2,ii*2)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*2,1)=global_stiffness_matrix(ii*2,ii*2)*displacement_left(2);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    % 右边界
    if boundary_condition_right(1)==1    % 如果右边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right    % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*2-1,ii*2-1)=global_stiffness_matrix(ii*2-1,ii*2-1)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*2-1,1)=global_stiffness_matrix(ii*2-1,ii*2-1)*displacement_right(1);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_right(2)==1    % 如果右边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right    % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*2,ii*2)=global_stiffness_matrix(ii*2,ii*2)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*2,1)=global_stiffness_matrix(ii*2,ii*2)*displacement_right(2);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    % 下边界
    if boundary_condition_bottom(1)==1    % 如果下边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_bottom    % 对下边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*2-1,ii*2-1)=global_stiffness_matrix(ii*2-1,ii*2-1)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*2-1,1)=global_stiffness_matrix(ii*2-1,ii*2-1)*displacement_bottom(1);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_bottom(2)==1    % 如果下边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_bottom    % 对下边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*2,ii*2)=global_stiffness_matrix(ii*2,ii*2)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*2,1)=global_stiffness_matrix(ii*2,ii*2)*displacement_bottom(2);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    % 上边界
    if boundary_condition_top(1)==1    % 如果上边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_top    % 对上边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*2-1,ii*2-1)=global_stiffness_matrix(ii*2-1,ii*2-1)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*2-1,1)=global_stiffness_matrix(ii*2-1,ii*2-1)*displacement_top(1);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_top(2)==1    % 如果上边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_top    % 对上边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*2,ii*2)=global_stiffness_matrix(ii*2,ii*2)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*2,1)=global_stiffness_matrix(ii*2,ii*2)*displacement_top(2);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加均布荷载条件
    % 左边界
    if boundary_condition_left(1)==2    % 如果左边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left    % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);    % ii表示左边界施加荷载条件的粒子编号
            global_load_vector(ii*2-1,1)=global_load_vector(ii*2-1,1)+uniform_load_left(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_left(2)==2    % 如果左边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left    % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);    % ii表示左边界施加荷载条件的粒子编号
            global_load_vector(ii*2,1)=global_load_vector(ii*2,1)+uniform_load_left(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    % 右边界
    if boundary_condition_right(1)==2    % 如果右边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right    % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);    % ii表示右边界施加荷载条件的粒子编号
            global_load_vector(ii*2-1,1)=global_load_vector(ii*2-1,1)+uniform_load_right(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_right(2)==2    % 如果右边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right    % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);    % ii表示右边界施加荷载条件的粒子编号
            global_load_vector(ii*2,1)=global_load_vector(ii*2,1)+uniform_load_right(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end 
    end
    % 下边界
    if boundary_condition_bottom(1)==2    % 如果下边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_bottom    % 对下边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_bottom(i);    % ii表示下边界施加荷载条件的粒子编号
            global_load_vector(ii*2-1,1)=global_load_vector(ii*2-1,1)+uniform_load_bottom(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_bottom(2)==2    % 如果下边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_bottom    % 对下边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_bottom(i);    % ii表示下边界施加荷载条件的粒子编号
            global_load_vector(ii*2,1)=global_load_vector(ii*2,1)+uniform_load_bottom(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    % 上边界
    if boundary_condition_top(1)==2    % 如果上边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_top    % 对上边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_top(i);    % ii表示上边界施加荷载条件的粒子编号
            global_load_vector(ii*2-1,1)=global_load_vector(ii*2-1,1)+uniform_load_top(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_top(2)==2    % 如果上边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_top    % 对上边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_top(i);    % ii表示上边界施加荷载条件的粒子编号
            global_load_vector(ii*2,1)=global_load_vector(ii*2,1)+uniform_load_top(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加体力
    if body_force(1)~=0    % 如果x方向体力密度不为0，说明施加了x方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle    % 对所有粒子循环
            global_load_vector(ii*2-1,1)=global_load_vector(ii*2-1,1)+body_force(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if body_force(2)~=0    % 如果y方向体力密度不为0，说明施加了x方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle    % 对所有粒子循环
            global_load_vector(ii*2,1)=global_load_vector(ii*2,1)+body_force(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%三维情况
if dimension==3
%%%%%%%%%%%%%%%%%%%%%%%
%施加位移边界条件
    % 左边界
    if boundary_condition_left(1)==1    % 如果左边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left    % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-2,ii*3-2)=global_stiffness_matrix(ii*3-2,ii*3-2)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-2,1)=global_stiffness_matrix(ii*3-2,ii*3-2)*displacement_left(1);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_left(2)==1    % 如果左边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left    % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-1,ii*3-1)=global_stiffness_matrix(ii*3-1,ii*3-1)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-1,1)=global_stiffness_matrix(ii*3-1,ii*3-1)*displacement_left(2);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_left(3)==1    % 如果左边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left    % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3,ii*3)=global_stiffness_matrix(ii*3,ii*3)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3,1)=global_stiffness_matrix(ii*3,ii*3)*displacement_left(3);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    % 右边界
    if boundary_condition_right(1)==1    % 如果右边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right    % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-2,ii*3-2)=global_stiffness_matrix(ii*3-2,ii*3-2)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-2,1)=global_stiffness_matrix(ii*3-2,ii*3-2)*displacement_right(1);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_right(2)==1    % 如果右边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right     % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-1,ii*3-1)=global_stiffness_matrix(ii*3-1,ii*3-1)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-1,1)=global_stiffness_matrix(ii*3-1,ii*3-1)*displacement_right(2);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_right(3)==1    % 如果右边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right    % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3,ii*3)=global_stiffness_matrix(ii*3,ii*3)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3,1)=global_stiffness_matrix(ii*3,ii*3)*displacement_right(3);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    % 下边界
    if boundary_condition_bottom(1)==1    % 如果下边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_bottom    % 对下边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-2,ii*3-2)=global_stiffness_matrix(ii*3-2,ii*3-2)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-2,1)=global_stiffness_matrix(ii*3-2,ii*3-2)*displacement_bottom(1);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_bottom(2)==1    % 如果下边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_bottom    % 对下边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-1,ii*3-1)=global_stiffness_matrix(ii*3-1,ii*3-1)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-1,1)=global_stiffness_matrix(ii*3-1,ii*3-1)*displacement_bottom(2);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_bottom(3)==1    % 如果下边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_bottom    % 对下边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3,ii*3)=global_stiffness_matrix(ii*3,ii*3)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3,1)=global_stiffness_matrix(ii*3,ii*3)*displacement_bottom(3);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    % 上边界
    if boundary_condition_top(1)==1    % 如果上边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_top    % 对上边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-2,ii*3-2)=global_stiffness_matrix(ii*3-2,ii*3-2)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-2,1)=global_stiffness_matrix(ii*3-2,ii*3-2)*displacement_top(1);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_top(2)==1    % 如果上边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_top    % 对上边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-1,ii*3-1)=global_stiffness_matrix(ii*3-1,ii*3-1)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-1,1)=global_stiffness_matrix(ii*3-1,ii*3-1)*displacement_top(2);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_top(3)==1    % 如果上边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_top    % 对上边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3,ii*3)=global_stiffness_matrix(ii*3,ii*3)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3,1)=global_stiffness_matrix(ii*3,ii*3)*displacement_top(3);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    % 后边界
    if boundary_condition_behind(1)==1    % 如果后边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_behind    % 对后边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_behind(i);    % ii表示后边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-2,ii*3-2)=global_stiffness_matrix(ii*3-2,ii*3-2)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-2,1)=global_stiffness_matrix(ii*3-2,ii*3-2)*displacement_behind(1);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_behind(2)==1    % 如果后边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_behind    % 对后边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_behind(i);    % ii表示后边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-1,ii*3-1)=global_stiffness_matrix(ii*3-1,ii*3-1)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-1,1)=global_stiffness_matrix(ii*3-1,ii*3-1)*displacement_behind(2);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_behind(3)==1    % 如果后边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_behind    % 对后边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_behind(i);    % ii表示后边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3,ii*3)=global_stiffness_matrix(ii*3,ii*3)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3,1)=global_stiffness_matrix(ii*3,ii*3)*displacement_behind(3);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    % 前边界
    if boundary_condition_front(1)==1    % 如果前边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_front    % 对前边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_front(i);    % ii表示前边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-2,ii*3-2)=global_stiffness_matrix(ii*3-2,ii*3-2)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-2,1)=global_stiffness_matrix(ii*3-2,ii*3-2)*displacement_front(1);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_front(2)==1    % 如果前边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_front    % 对前边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_front(i);    % ii表示前边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3-1,ii*3-1)=global_stiffness_matrix(ii*3-1,ii*3-1)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3-1,1)=global_stiffness_matrix(ii*3-1,ii*3-1)*displacement_front(2);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
    if boundary_condition_front(3)==1    % 如果前边界定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_front    % 对前边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_front(i);    % ii表示前边界施加位移边界条件的粒子编号
            global_stiffness_matrix(ii*3,ii*3)=global_stiffness_matrix(ii*3,ii*3)*1.0e8;    % 采用乘大数法对总体刚度矩阵进行操作
            global_load_vector(ii*3,1)=global_stiffness_matrix(ii*3,ii*3)*displacement_front(3);    % 采用乘大数法对总体荷载向量进行操作
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加均布荷载条件
    % 左边界
    if boundary_condition_left(1)==2    % 如果左边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left    % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-2,1)=global_load_vector(ii*3-2,1)+uniform_load_left(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_left(2)==2    % 如果左边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left    % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-1,1)=global_load_vector(ii*3-1,1)+uniform_load_left(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_left(3)==2    % 如果左边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left    % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_load_vector(ii*3,1)=global_load_vector(ii*3,1)+uniform_load_left(3)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    % 右边界
    if boundary_condition_right(1)==2    % 如果右边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right    % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-2,1)=global_load_vector(ii*3-2,1)+uniform_load_right(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_right(2)==2    % 如果右边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right    % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-1,1)=global_load_vector(ii*3-1,1)+uniform_load_right(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_right(3)==2    % 如果右边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right    % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_load_vector(ii*3,1)=global_load_vector(ii*3,1)+uniform_load_right(3)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    % 下边界
    if boundary_condition_bottom(1)==2    % 如果下边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_bottom    % 对下边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-2,1)=global_load_vector(ii*3-2,1)+uniform_load_bottom(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_bottom(2)==2    % 如果下边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_bottom     % 对下边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-1,1)=global_load_vector(ii*3-1,1)+uniform_load_bottom(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_bottom(3)==2    % 如果下边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_bottom    % 对下边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_load_vector(ii*3,1)=global_load_vector(ii*3,1)+uniform_load_bottom(3)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    % 上边界
    if boundary_condition_top(1)==2    % 如果上边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_top    % 对上边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-2,1)=global_load_vector(ii*3-2,1)+uniform_load_top(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_top(2)==2    % 如果上边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_top    % 对上边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-1,1)=global_load_vector(ii*3-1,1)+uniform_load_top(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_top(3)==2    % 如果上边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_top     % 对上边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_vector(ii*3,1)=global_load_vector(ii*3,1)+uniform_load_top(3)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    % 后边界
    if boundary_condition_behind(1)==2    % 如果后边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_behind    % 对后边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_behind(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-2,1)=global_load_vector(ii*3-2,1)+uniform_load_behind(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_behind(2)==2    % 如果后边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_behind    % 对后边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_behind(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-1,1)=global_load_vector(ii*3-1,1)+uniform_load_behind(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end  
    end
    if boundary_condition_behind(3)==2    % 如果后边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_behind    % 对后边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_behind(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_vector(ii*3,1)=global_load_vector(ii*3,1)+uniform_load_behind(3)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    % 前边界
    if boundary_condition_front(1)==2    % 如果前边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_front    % 对前边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_front(i);    % ii表示前边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-2,1)=global_load_vector(ii*3-2,1)+uniform_load_front(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_front(2)==2    % 如果前边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_front    % 对前边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_front(i);    % ii表示前边界施加位移边界条件的粒子编号
            global_load_vector(ii*3-1,1)=global_load_vector(ii*3-1,1)+uniform_load_front(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if boundary_condition_front(3)==2    % 如果前边界定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_front    % 对前边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_front(i);    % ii表示前边界施加位移边界条件的粒子编号
            global_load_vector(ii*3,1)=global_load_vector(ii*3,1)+uniform_load_front(3)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加体力
    if body_force(1)~=0    % 如果x方向体力密度不为0，说明施加了x方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle     % 对所有粒子循环
            global_load_vector(ii*3-2,1)=global_load_vector(ii*3-2,1)+body_force(1)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if body_force(2)~=0    % 如果y方向体力密度不为0，说明施加了x方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle    % 对所有粒子循环
            global_load_vector(ii*3-1,1)=global_load_vector(ii*3-1,1)+body_force(2)*volume_particle;    % 在总体荷载向量对应位置增加荷载
        end
    end
    if body_force(3)~=0    % 如果z方向体力密度不为0，说明施加了x方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle    % 对所有粒子循环
            global_load_vector(ii*3,1)=global_load_vector(ii*3,1)+body_force(3)*volume_particle;     % 在总体荷载向量对应位置增加荷载
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global_displacement_vector=global_stiffness_matrix\global_load_vector;  % 求解总体位移荷载
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%