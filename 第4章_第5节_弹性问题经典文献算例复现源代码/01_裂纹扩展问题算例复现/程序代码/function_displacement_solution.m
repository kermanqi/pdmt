function[global_velocity_vector,global_displacement_vector]=function_displacement_solution(dimension,size_particle,size_displacement_boundary_left,particle_boundary_left, ...
         size_displacement_boundary_right,particle_boundary_right,size_displacement_boundary_bottom,particle_boundary_bottom,size_displacement_boundary_top, ...
         particle_boundary_top,size_displacement_boundary_behind,particle_boundary_behind,size_displacement_boundary_front,particle_boundary_front, ...
         size_uniform_load_left,particle_uniform_load_left,size_uniform_load_right,particle_uniform_load_right,size_uniform_load_bottom,particle_uniform_load_bottom, ...
         size_uniform_load_top,particle_uniform_load_top,size_uniform_load_behind,particle_uniform_load_behind,size_uniform_load_front,particle_uniform_load_front, ...
         velocity_left,velocity_right,velocity_bottom,velocity_top,velocity_behind,velocity_front,uniform_load_left,uniform_load_right, ...
         uniform_load_bottom,uniform_load_top,uniform_load_behind,uniform_load_front,body_force,boundary_condition_left,boundary_condition_right,boundary_condition_bottom, ...
         boundary_condition_top,boundary_condition_behind,boundary_condition_front,itime,dt,output_interval,particle,horizon_particle,size_horizon_particle,c0,density, ...
         global_velocity_vector,global_displacement_vector,volume_particle,damage_particle)


itime=double(itime);  % 将计算步数改写为浮点型变量
global_load_density_vector=zeros(size_particle*dimension,1); % 总体荷载密度矢量赋初值
global_accelerate_vector=zeros(size_particle*dimension,1);  % 总体加速度矢量赋初值
global_velocity0_vector=global_velocity_vector;  % 将上一时间步计算得到的总体速度矢量赋值给 global_velocity0_vector
global_displacement0_vector=global_displacement_vector;  % 将上一时间步计算得到的总体位移矢量赋值给 global_displacement0_vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%一维情况
if dimension==1
%%%%%%%%%%%%%%%%%%%%%%%
%施加位移边界条件
    % 左边界
    if boundary_condition_left(1)==1                                % 如果左边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left                     % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);                           % ii表示左边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii,1)=velocity_left(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii,1)=velocity_left(1)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    % 右边界
    if boundary_condition_right(1)==1                           % 如果右边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right                % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);                      % ii表示右边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii,1)=velocity_right(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii,1)=velocity_right(1)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加均布荷载条件
    % 左边界
    if boundary_condition_left(1)==2    % 如果左边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left  % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);   % ii表示左边界施加荷载条件的粒子编号
            global_load_density_vector(ii,1)=global_load_density_vector(ii,1)+uniform_load_left(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    % 右边界
    if boundary_condition_right(1)==2    % 如果右边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right  % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);   % ii表示右边界施加荷载条件的粒子编号
            global_load_density_vector(ii,1)=global_load_density_vector(ii,1)+uniform_load_right(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加体力
    if body_force(1)~=0    % 如果x方向体力密度不为0，说明施加了x方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle % 对所有粒子循环
            global_load_density_vector(ii,1)=global_load_density_vector(ii,1)+body_force(1); % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%位移求解
    for i=1:size_particle    % 对所有粒子循环
        force_density_x=0;    % x方向力密度矢量清零
        xx_i=particle(i,1);     % 获取第i个粒子的x方向坐标
        uu_i=global_displacement0_vector(i*2-1,1);     % 获取当前时刻第i个粒子的x方向位移
        for j=1:size_horizon_particle       % 对第i个粒子作用范围内的所有粒子进行循环
            jj=horizon_particle(i,j);       % 第i个粒子作用范围内的第j个粒子为粒子jj
            if jj>0       % 如果第i个粒子作用范围内的第j个粒子存在，则执行循环条件内的语句
                if damage_particle(i,j)<0.5     % 如果第i个粒子与第jj个粒子组成的键没有发生破坏，则执行条件语句中的内容
                    xx_j=particle(jj,1);     % 获取第jj个粒子的x方向坐标
                    uu_j=global_displacement0_vector(jj*2-1,1);     % 获取当前时刻第jj个粒子的x方向位移
                    length1_i_jj=abs(xx_j-xx_i);   % 初始构型下粒子i与粒子jj之间的距离为length_i_jj
                    length2_i_jj=abs(xx_j+uu_j-xx_i-uu_i);   % 当前构型下粒子i与粒子jj之间的距离为length2_i_jj
                    s0=(length2_i_jj-length1_i_jj)/length1_i_jj;   % 粒子i与粒子jj之间键的伸长率为s0
                    force_density_x=force_density_x+c0*s0*(xx_j-xx_i)/length1_i_jj; % 粒子i与粒子jj在x方向的作用加到粒子i所受到的x方向力密度矢量中
                end
            end
        end
        global_accelerate_vector(i,1)=(force_density_x+global_load_density_vector(i,1)*volume_particle)/density; % 当前时刻粒子i在x方向的加速度
    end
    global_velocity_vector=global_velocity0_vector+global_accelerate_vector*dt;          % 更新当前时刻所有粒子的速度
    global_displacement_vector=global_displacement0_vector+global_velocity_vector*dt;          % 更新当前时刻所有粒子的位移
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%二维情况
if dimension==2
%%%%%%%%%%%%%%%%%%%%%%%
%施加位移边界条件
    % 左边界
    if boundary_condition_left(1)==1    % 如果左边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left    % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*2-1,1)=velocity_left(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*2-1,1)=velocity_left(1)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_left(2)==1    % 如果左边界y方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left    % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*2,1)=velocity_left(2);     % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*2,1)=velocity_left(2)*itime*dt;      % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    % 右边界
    if boundary_condition_right(1)==1    % 如果右边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right    % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*2-1,1)=velocity_right(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*2-1,1)=velocity_right(1)*itime*dt;      % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_right(2)==1    % 如果右边界y方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right    % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*2,1)=velocity_right(2);     % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*2,1)=velocity_right(2)*itime*dt;      % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    % 下边界
    if boundary_condition_bottom(1)==1    % 如果下边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_bottom    % 对下边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*2-1,1)=velocity_bottom(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*2-1,1)=velocity_bottom(1)*itime*dt;      % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_bottom(2)==1    % 如果下边界y方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_bottom    % 对下边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*2,1)=velocity_bottom(2);     % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*2,1)=velocity_bottom(2)*itime*dt;      % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    % 上边界
    if boundary_condition_top(1)==1    % 如果上边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_top    % 对上边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*2-1,1)=velocity_top(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*2-1,1)=velocity_top(1)*itime*dt;      % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_top(2)==1    % 如果上边界y方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_top    % 对上边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*2,1)=velocity_top(2);     % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*2,1)=velocity_top(2)*itime*dt;      % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加均布荷载条件
    % 左边界
    if boundary_condition_left(1)==2    % 如果左边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left    % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);    % ii表示左边界施加荷载条件的粒子编号
            global_load_density_vector(ii*2-1,1)=global_load_density_vector(ii*2-1,1)+uniform_load_left(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_left(2)==2    % 如果左边界y方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left    % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);    % ii表示左边界施加荷载条件的粒子编号
            global_load_density_vector(ii*2,1)=global_load_density_vector(ii*2,1)+uniform_load_left(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    % 右边界
    if boundary_condition_right(1)==2    % 如果右边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right    % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);    % ii表示右边界施加荷载条件的粒子编号
            global_load_density_vector(ii*2-1,1)=global_load_density_vector(ii*2-1,1)+uniform_load_right(1);    % 在总体荷密度载向量对应位置增加荷载密度
        end
    end
    if boundary_condition_right(2)==2    % 如果右边界y方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right    % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);    % ii表示右边界施加荷载条件的粒子编号
            global_load_density_vector(ii*2,1)=global_load_density_vector(ii*2,1)+uniform_load_right(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end 
    end
    % 下边界
    if boundary_condition_bottom(1)==2    % 如果下边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_bottom    % 对下边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_bottom(i);    % ii表示下边界施加荷载条件的粒子编号
            global_load_density_vector(ii*2-1,1)=global_load_density_vector(ii*2-1,1)+uniform_load_bottom(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_bottom(2)==2    % 如果下边界y方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_bottom    % 对下边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_bottom(i);    % ii表示下边界施加荷载条件的粒子编号
            global_load_density_vector(ii*2,1)=global_load_density_vector(ii*2,1)+uniform_load_bottom(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    % 上边界
    if boundary_condition_top(1)==2    % 如果上边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_top    % 对上边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_top(i);    % ii表示上边界施加荷载条件的粒子编号
            global_load_density_vector(ii*2-1,1)=global_load_density_vector(ii*2-1,1)+uniform_load_top(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_top(2)==2    % 如果上边界y方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_top    % 对上边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_top(i);    % ii表示上边界施加荷载条件的粒子编号
            global_load_density_vector(ii*2,1)=global_load_density_vector(ii*2,1)+uniform_load_top(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加体力
    if body_force(1)~=0    % 如果x方向体力密度不为0，说明施加了x方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle    % 对所有粒子循环
            global_load_density_vector(ii*2-1,1)=global_load_density_vector(ii*2-1,1)+body_force(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if body_force(2)~=0    % 如果y方向体力密度不为0，说明施加了y方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle    % 对所有粒子循环
            global_load_density_vector(ii*2,1)=global_load_density_vector(ii*2,1)+body_force(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%位移求解
    for i=1:size_particle    % 对所有粒子循环
        force_density_x=0;    % x方向力密度矢量清零
        force_density_y=0;    % y方向力密度矢量清零
        xx_i=particle(i,1);     % 获取第i个粒子的x方向坐标
        yy_i=particle(i,2);     % 获取第i个粒子的y方向坐标
        uu_i=global_displacement0_vector(i*2-1,1);     % 获取当前时刻第i个粒子的x方向位移
        vv_i=global_displacement0_vector(i*2,1);       % 获取当前时刻第i个粒子的y方向位移
        for j=1:size_horizon_particle       % 对第i个粒子作用范围内的所有粒子进行循环
            jj=horizon_particle(i,j);       % 第i个粒子作用范围内的第j个粒子为粒子jj
            if jj>0       % 如果第i个粒子作用范围内的第j个粒子存在，则执行循环条件内的语句
                if damage_particle(i,j)<0.5     % 如果第i个粒子与第jj个粒子组成的键没有发生破坏，则执行条件语句中的内容
                    xx_j=particle(jj,1);     % 获取第jj个粒子的x方向坐标
                    yy_j=particle(jj,2);     % 获取第jj个粒子的y方向坐标
                    uu_j=global_displacement0_vector(jj*2-1,1);     % 获取当前时刻第jj个粒子的x方向位移
                    vv_j=global_displacement0_vector(jj*2,1);     % 获取当前时刻第jj个粒子的y方向位移
                    length1_i_jj=sqrt((xx_j-xx_i)*(xx_j-xx_i)+(yy_j-yy_i)*(yy_j-yy_i));   % 初始构型下粒子i与粒子jj之间的距离为length_i_jj
                    length2_i_jj=sqrt((xx_j+uu_j-xx_i-uu_i)*(xx_j+uu_j-xx_i-uu_i)+(yy_j+vv_j-yy_i-vv_i)*(yy_j+vv_j-yy_i-vv_i));   % 当前构型下粒子i与粒子jj之间的距离为length2_i_jj
                    s0=(length2_i_jj-length1_i_jj)/length1_i_jj;   % 粒子i与粒子jj之间键的伸长率为s0
                    force_density_x=force_density_x+c0*s0*(xx_j-xx_i)/length1_i_jj; % 粒子i与粒子jj在x方向的作用加到粒子i所受到的x方向力密度矢量中
                    force_density_y=force_density_y+c0*s0*(yy_j-yy_i)/length1_i_jj; % 粒子i与粒子jj在y方向的作用加到粒子i所受到的y方向力密度矢量中
                end
            end
        end
        global_accelerate_vector(i*2-1,1)=(force_density_x+global_load_density_vector(i*2-1,1)*volume_particle)/density; % 当前时刻粒子i在x方向的加速度
        global_accelerate_vector(i*2,1)=(force_density_y+global_load_density_vector(i*2,1)*volume_particle)/density; % 当前时刻粒子i在y方向的加速度
    end
    global_velocity_vector=global_velocity0_vector+global_accelerate_vector*dt;          % 更新当前时刻所有粒子的速度
    global_displacement_vector=global_displacement0_vector+global_velocity_vector*dt;          % 更新当前时刻所有粒子的位移
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%三维情况
if dimension==3
%%%%%%%%%%%%%%%%%%%%%%%
%施加位移边界条件
    % 左边界
    if boundary_condition_left(1)==1    % 如果左边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left    % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-2,1)=velocity_left(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-2,1)=velocity_left(1)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_left(2)==1    % 如果左边界y方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left    % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-1,1)=velocity_left(2);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-1,1)=velocity_left(2)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_left(3)==1    % 如果左边界z方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left    % 对左边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3,1)=velocity_left(3);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3,1)=velocity_left(3)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    % 右边界
    if boundary_condition_right(1)==1    % 如果右边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right    % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-2,1)=velocity_right(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-2,1)=velocity_right(1)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_right(2)==1    % 如果右边界y方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right     % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-1,1)=velocity_right(2);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-1,1)=velocity_right(2)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_right(3)==1    % 如果右边界z方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right    % 对右边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_right(i);    % ii表示右边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3,1)=velocity_right(3);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3,1)=velocity_right(3)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    % 下边界
    if boundary_condition_bottom(1)==1    % 如果下边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_bottom    % 对下边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-2,1)=velocity_bottom(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-2,1)=velocity_bottom(1)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_bottom(2)==1    % 如果下边界y方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_bottom    % 对下边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-1,1)=velocity_bottom(2);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-1,1)=velocity_bottom(2)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_bottom(3)==1    % 如果下边界z方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_bottom    % 对下边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3,1)=velocity_bottom(3);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3,1)=velocity_bottom(3)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    % 上边界
    if boundary_condition_top(1)==1    % 如果上边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_top    % 对上边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-2,1)=velocity_top(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-2,1)=velocity_top(1)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_top(2)==1    % 如果上边界y方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_top    % 对上边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-1,1)=velocity_top(2);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-1,1)=velocity_top(2)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_top(3)==1    % 如果上边界z方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_top    % 对上边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3,1)=velocity_top(3);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3,1)=velocity_top(3)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    % 后边界
    if boundary_condition_behind(1)==1    % 如果后边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_behind    % 对后边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_behind(i);    % ii表示后边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-2,1)=velocity_behind(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-2,1)=velocity_behind(1)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_behind(2)==1    % 如果后边界y方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_behind    % 对后边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_behind(i);    % ii表示后边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-1,1)=velocity_behind(2);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-1,1)=velocity_behind(2)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_behind(3)==1    % 如果后边界z方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_behind    % 对后边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_behind(i);    % ii表示后边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3,1)=velocity_behind(3);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3,1)=velocity_behind(3)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    % 前边界
    if boundary_condition_front(1)==1    % 如果前边界x方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_front    % 对前边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_front(i);    % ii表示前边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-2,1)=velocity_front(1);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-2,1)=velocity_front(1)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_front(2)==1    % 如果前边界y方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_front    % 对前边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_front(i);    % ii表示前边界施加位移边界条件的粒子编号
            global_velocity0_vector(ii*3-1,1)=velocity_front(2);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3-1,1)=velocity_front(2)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
    if boundary_condition_front(3)==1    % 如果前边界z方向定义了位移边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_front    % 对前边界施加位移边界条件的所有粒子循环
            ii=particle_boundary_front(i);    % ii表示前边界施加位移边界条件的粒子编号
           global_velocity0_vector(ii*3,1)=velocity_front(3);    % 在总体速度向量的对应位置赋予速度边界条件
            global_displacement0_vector(ii*3,1)=velocity_front(3)*itime*dt;     % 在总体位移向量的对应位置赋予该时间步的位移
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加均布荷载条件
    % 左边界
    if boundary_condition_left(1)==2    % 如果左边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left    % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3-2,1)=global_load_density_vector(ii*3-2,1)+uniform_load_left(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_left(2)==2    % 如果左边界y方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left    % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3-1,1)=global_load_density_vector(ii*3-1,1)+uniform_load_left(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_left(3)==2    % 如果左边界z方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_left    % 对左边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_left(i);    % ii表示左边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3,1)=global_load_density_vector(ii*3,1)+uniform_load_left(3);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    % 右边界
    if boundary_condition_right(1)==2    % 如果右边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right    % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);    % ii表示右边界施加位移边界条件的粒子编号
             global_load_density_vector(ii*3-2,1)=global_load_density_vector(ii*3-2,1)+uniform_load_right(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_right(2)==2    % 如果右边界y方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right    % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);    % ii表示右边界施加位移边界条件的粒子编号
             global_load_density_vector(ii*3-1,1)=global_load_density_vector(ii*3-1,1)+uniform_load_right(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_right(3)==2    % 如果右边界z方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_right    % 对右边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_right(i);    % ii表示右边界施加位移边界条件的粒子编号
             global_load_density_vector(ii*3,1)=global_load_density_vector(ii*3,1)+uniform_load_right(3);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    % 下边界
    if boundary_condition_bottom(1)==2    % 如果下边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_bottom    % 对下边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3-2,1)=global_load_density_vector(ii*3-2,1)+uniform_load_bottom(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_bottom(2)==2    % 如果下边界y方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_bottom     % 对下边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3-1,1)=global_load_density_vector(ii*3-1,1)+uniform_load_bottom(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_bottom(3)==2    % 如果下边界z方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_bottom    % 对下边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_bottom(i);    % ii表示下边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3,1)=global_load_density_vector(ii*3,1)+uniform_load_bottom(3);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    % 上边界
    if boundary_condition_top(1)==2    % 如果上边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_top    % 对上边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3-2,1)=global_load_density_vector(ii*3-2,1)+uniform_load_top(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_top(2)==2    % 如果上边界y方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_top    % 对上边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3-1,1)=global_load_density_vector(ii*3-1,1)+uniform_load_top(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_top(3)==2    % 如果上边界z方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_top     % 对上边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_top(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3,1)=global_load_density_vector(ii*3,1)+uniform_load_top(3);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    % 后边界
    if boundary_condition_behind(1)==2    % 如果后边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_behind    % 对后边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_behind(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3-2,1)=global_load_density_vector(ii*3-2,1)+uniform_load_behind(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_behind(2)==2    % 如果后边界y方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_behind    % 对后边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_behind(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3-1,1)=global_load_density_vector(ii*3-1,1)+uniform_load_behind(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end  
    end
    if boundary_condition_behind(3)==2    % 如果后边界z方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_behind    % 对后边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_behind(i);    % ii表示上边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3,1)=global_load_density_vector(ii*3,1)+uniform_load_behind(3);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    % 前边界
    if boundary_condition_front(1)==2    % 如果前边界x方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_front    % 对前边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_front(i);    % ii表示前边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3-2,1)=global_load_density_vector(ii*3-2,1)+uniform_load_front(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_front(2)==2    % 如果前边界y方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_front    % 对前边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_front(i);    % ii表示前边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3-1,1)=global_load_density_vector(ii*3-1,1)+uniform_load_front(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if boundary_condition_front(3)==2    % 如果前边界z方向定义了荷载条件，则执行条件语句中的内容
        for i=1:size_uniform_load_front    % 对前边界施加荷载条件的所有粒子循环
            ii=particle_uniform_load_front(i);    % ii表示前边界施加位移边界条件的粒子编号
            global_load_density_vector(ii*3,1)=global_load_density_vector(ii*3,1)+uniform_load_front(3);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%施加体力
    if body_force(1)~=0    % 如果x方向体力密度不为0，说明施加了x方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle     % 对所有粒子循环
            global_load_density_vector(ii*3-2,1)=global_load_density_vector(ii*3-2,1)+body_force(1);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if body_force(2)~=0    % 如果y方向体力密度不为0，说明施加了y方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle    % 对所有粒子循环
            global_load_density_vector(ii*3-1,1)=global_load_density_vector(ii*3-1,1)+body_force(2);    % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
    if body_force(3)~=0    % 如果z方向体力密度不为0，说明施加了x方向的体力荷载，则执行条件语句中的内容
        for i=1:size_particle    % 对所有粒子循环
            global_load_density_vector(ii*3,1)=global_load_density_vector(ii*3,1)+body_force(3);     % 在总体荷载密度向量对应位置增加荷载密度
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
%位移求解
    for i=1:size_particle    % 对所有粒子循环
        force_density_x=0;    % x方向力密度矢量清零
        force_density_y=0;    % y方向力密度矢量清零
        force_density_z=0;    % z方向力密度矢量清零
        xx_i=particle(i,1);     % 获取第i个粒子的x方向坐标
        yy_i=particle(i,2);     % 获取第i个粒子的y方向坐标
        zz_i=particle(i,3);     % 获取第i个粒子的z方向坐标
        uu_i=global_displacement0_vector(i*3-2,1);     % 获取当前时刻第i个粒子的x方向位移
        vv_i=global_displacement0_vector(i*3-1,1);       % 获取当前时刻第i个粒子的y方向位移
        ww_i=global_displacement0_vector(i*3,1);       % 获取当前时刻第i个粒子的z方向位移
        for j=1:size_horizon_particle       % 对第i个粒子作用范围内的所有粒子进行循环
            jj=horizon_particle(i,j);       % 第i个粒子作用范围内的第j个粒子为粒子jj
            if jj>0       % 如果第i个粒子作用范围内的第j个粒子存在，则执行循环条件内的语句
                if damage_particle(i,j)<0.5     % 如果第i个粒子与第jj个粒子组成的键没有发生破坏，则执行条件语句中的内容
                    xx_j=particle(jj,1);     % 获取第jj个粒子的x方向坐标
                    yy_j=particle(jj,2);     % 获取第jj个粒子的y方向坐标
                    zz_j=particle(jj,3);     % 获取第jj个粒子的z方向坐标
                    uu_j=global_displacement0_vector(jj*3-2,1);     % 获取当前时刻第jj个粒子的x方向位移
                    vv_j=global_displacement0_vector(jj*3-1,1);     % 获取当前时刻第jj个粒子的y方向位移
                    ww_j=global_displacement0_vector(jj*3,1);     % 获取当前时刻第jj个粒子的z方向位移
                    length1_i_jj=sqrt((xx_j-xx_i)*(xx_j-xx_i)+(yy_j-yy_i)*(yy_j-yy_i)+(zz_j-zz_i)*(zz_j-zz_i));   % 初始构型下粒子i与粒子jj之间的距离为length_i_jj
                    length2_i_jj=sqrt((xx_j+uu_j-xx_i-uu_i)*(xx_j+uu_j-xx_i-uu_i)+(yy_j+vv_j-yy_i-vv_i)*(yy_j+vv_j-yy_i-vv_i)+(zz_j+ww_j-zz_i-ww_i)*(zz_j+ww_j-zz_i-ww_i));   % 当前构型下粒子i与粒子jj之间的距离为length2_i_jj
                    s0=(length2_i_jj-length1_i_jj)/length1_i_jj;   % 粒子i与粒子jj之间键的伸长率为s0
                    force_density_x=force_density_x+c0*s0*(xx_j-xx_i)/length1_i_jj; % 粒子i与粒子jj在x方向的作用加到粒子i所受到的x方向力密度矢量中
                    force_density_y=force_density_y+c0*s0*(yy_j-yy_i)/length1_i_jj; % 粒子i与粒子jj在y方向的作用加到粒子i所受到的y方向力密度矢量中
                    force_density_z=force_density_z+c0*s0*(zz_j-zz_i)/length1_i_jj; % 粒子i与粒子jj在z方向的作用加到粒子i所受到的z方向力密度矢量中
                end
            end
        end
        global_accelerate_vector(i*3-2,1)=(force_density_x+global_load_density_vector(i*3-2,1)*volume_particle)/density; % 当前时刻粒子i在x方向的加速度
        global_accelerate_vector(i*3-1,1)=(force_density_y+global_load_density_vector(i*3-1,1)*volume_particle)/density; % 当前时刻粒子i在y方向的加速度
        global_accelerate_vector(i*3,1)=(force_density_y+global_load_density_vector(i*3,1)*volume_particle)/density; % 当前时刻粒子i在z方向的加速度
    end
    global_velocity_vector=global_velocity0_vector+global_accelerate_vector*dt;          % 更新当前时刻所有粒子的速度
    global_displacement_vector=global_displacement0_vector+global_velocity_vector*dt;          % 更新当前时刻所有粒子的位移
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(itime,output_interval)==0  % 如果itime对output_interval取余数为0，则执行条件语句中的内容
    str1='displacement_';           % 定义输出文件名称中的第1部分 displacement_
    str2=num2str(itime);           % 定义输出文件名称中的第2部分 itime
    str3='.txt';                    % 定义输出文件名称中的第3部分 .txt，即文件拓展名
    str4=strcat(str1,str2,str3);    % 将文件名称写到一个变量中，即displacement_itime.txt，其中itime指具体计算步数，如果itime=1,则输出文件名为displacement_1.txt
    dlmwrite(str4,global_displacement_vector);  % 调用 dlmwrite 函数将global_displacement_vector写入文本文件displacement_itime.txt
end