function[global_velocity_vector,global_displacement_vector]=function_displacement_solution(dimension,size_particle,size_displacement_boundary_left,particle_boundary_left, ...
         size_displacement_boundary_right,particle_boundary_right,size_displacement_boundary_bottom,particle_boundary_bottom,size_displacement_boundary_top, ...
         particle_boundary_top, ...
         size_uniform_load_left,particle_uniform_load_left,size_uniform_load_right,particle_uniform_load_right,size_uniform_load_bottom,particle_uniform_load_bottom, ...
         size_uniform_load_top,particle_uniform_load_top, ...
         velocity_left,velocity_right,velocity_bottom,velocity_top,uniform_load_left,uniform_load_right, ...
         uniform_load_bottom,uniform_load_top,body_force,boundary_condition_left,boundary_condition_right,boundary_condition_bottom, ...
         boundary_condition_top,itime,dt,output_interval,particle,horizon_particle,size_horizon_particle,cf,cm,bn,bs,density, ...
         global_velocity_vector,global_displacement_vector,volume_particle,classify_bond,horizon_particle_n,size_horizon_particle_n, ...
         horizon_particle_s_lower,horizon_particle_s_upper,detla_n,detla_s,particle_ply)


itime=double(itime);  % 将计算步数改写为浮点型变量
global_load_density_vector=zeros(size_particle*dimension,1); % 总体荷载密度矢量赋初值
global_accelerate_vector=zeros(size_particle*dimension,1);  % 总体加速度矢量赋初值
global_velocity0_vector=global_velocity_vector;  % 将上一时间步计算得到的总体速度矢量赋值给 global_velocity0_vector
global_displacement0_vector=global_displacement_vector;  % 将上一时间步计算得到的总体位移矢量赋值给 global_displacement0_vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% 位移求解
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 面内作用
    for j=1:size_horizon_particle       % 对第i个粒子作用范围内的所有粒子进行循环
        jj=horizon_particle(i,j);       % 第i个粒子作用范围内的第j个粒子为粒子jj
        if jj>0       % 如果第i个粒子作用范围内的第j个粒子存在，则执行循环条件内的语句
        	xx_j=particle(jj,1);     % 获取第jj个粒子的x方向坐标
            yy_j=particle(jj,2);     % 获取第jj个粒子的y方向坐标
            uu_j=global_displacement0_vector(jj*3-2,1);     % 获取当前时刻第jj个粒子的x方向位移
            vv_j=global_displacement0_vector(jj*3-1,1);     % 获取当前时刻第jj个粒子的y方向位移
            length1_i_jj=sqrt((xx_j-xx_i)*(xx_j-xx_i)+(yy_j-yy_i)*(yy_j-yy_i));   % 初始构型下粒子i与粒子jj之间的距离为length_i_jj
            length2_i_jj=sqrt((xx_j+uu_j-xx_i-uu_i)*(xx_j+uu_j-xx_i-uu_i)+(yy_j+vv_j-yy_i-vv_i)*(yy_j+vv_j-yy_i-vv_i));   % 当前构型下粒子i与粒子jj之间的距离为length2_i_jj
            s0=(length2_i_jj-length1_i_jj)/length1_i_jj;   % 粒子i与粒子jj之间键的伸长率为s0
            if classify_bond(i,j)==1	% classify_bond(i,j)等于1表示纤维方向的键
                ii=particle_ply(i);     % 获取粒子i所在层数
                c0=cf(ii)+cm;   % 微模量系数取纤维键微模量系数与基体键微模量系数的和
            else	% classify_bond(i,j)不等于1表示其他方向的键
                c0=cm;   % 微模量系数取基体键微模量系数
            end
            force_density_x=force_density_x+c0*s0*(xx_j-xx_i)/length1_i_jj; % 粒子i与粒子jj在x方向的作用加到粒子i所受到的x方向力密度矢量中
            force_density_y=force_density_y+c0*s0*(yy_j-yy_i)/length1_i_jj; % 粒子i与粒子jj在y方向的作用加到粒子i所受到的y方向力密度矢量中
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 面外法向作用
	for j=1:size_horizon_particle_n           % 对所有与粒子i进行面外法向作用的粒子进行循环
        jj=horizon_particle_n(i,j);           % 与粒子i进行面外法向作用的第j个粒子的粒子编号为jj
        if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
            zz_j=particle(jj,3);    % 获取第jj个粒子的z坐标
            ww_j=global_displacement0_vector(jj*3,1);     % 获取当前时刻第jj个粒子的z方向位移
            length1_i_jj=abs(zz_j-zz_i);   % 初始构型下粒子i与粒子jj之间的距离为length1_i_jj
            length2_i_jj=abs(zz_j+ww_j-zz_i-ww_i);   % 当前构型下粒子i与粒子jj之间的距离为length2_i_jj
            s0=(length2_i_jj-length1_i_jj)/length1_i_jj;   % 粒子i与粒子jj之间键的伸长率为s0
            force_density_z=force_density_z+4*bn*detla_n*s0*(zz_j-zz_i)/length1_i_jj; % 粒子i与粒子jj在z方向的作用加到粒子i所受到的z方向力密度矢量中
        end
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 层间剪切作用
    for ii=1:size_horizon_particle_n    % 对所有与粒子i进行面外法向作用的粒子进行循环
        iin=horizon_particle_n(i,ii);   % 与粒子i进行面外法向作用的第ii个粒子的粒子编号为iin
        if iin>0                        % 如果第i个粒子域内的第ii个粒子存在，则执行条件语句中的内容
            xx_iin=particle(iin,1);  % 获取第iin个粒子的x坐标                   
            yy_iin=particle(iin,2);  % 获取第iin个粒子的y坐标                   
            zz_iin=particle(iin,3);  % 获取第iin个粒子的z坐标
            uu_iin=global_displacement0_vector(iin*3-2,1);     % 获取当前时刻第iin个粒子的x方向位移
            vv_iin=global_displacement0_vector(iin*3-1,1);       % 获取当前时刻第iin个粒子的y方向位移
            ww_iin=global_displacement0_vector(iin*3,1);       % 获取当前时刻第iin个粒子的z方向位移
            for j=1:size_horizon_particle           % 对所有与粒子i进行面内作用的粒子进行循环
                jj=horizon_particle(i,j);           % 与粒子i进行面内作用的第j个粒子的粒子编号为jj
                if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
                    xx_j=particle(jj,1);            % 获取第jj个粒子的x坐标
                    yy_j=particle(jj,2);            % 获取第jj个粒子的y坐标
                    zz_j=particle(jj,3);            % 获取第jj个粒子的z坐标
                    uu_j=global_displacement0_vector(jj*3-2,1);     % 获取当前时刻第jj个粒子的x方向位移
                    vv_j=global_displacement0_vector(jj*3-1,1);     % 获取当前时刻第jj个粒子的y方向位移
                    ww_j=global_displacement0_vector(jj*3,1);     % 获取当前时刻第jj个粒子的z方向位移
                    if ii==1                        % 如果ii等于1，则表示粒子iin为粒子i下面一层的粒子
                        jjs=horizon_particle_s_lower(i,j); % 获取位于粒子jj下面一层的粒子jjs
                    end
                    if ii==2                        % 如果ii等于2，则表示粒子iin为粒子i上面一层的粒子
                        jjs=horizon_particle_s_upper(i,j); % 获取位于粒子jj上面一层的粒子jjs
                    end
                    xx_jjs=particle(jjs,1);     % 获取第jjs个粒子的x坐标                
                    yy_jjs=particle(jjs,2);     % 获取第jjs个粒子的y坐标               
                    zz_jjs=particle(jjs,3);     % 获取第jjs个粒子的z坐标
                    uu_jjs=global_displacement0_vector(jjs*3-2,1);     % 获取当前时刻第jjs个粒子的x方向位移
                    vv_jjs=global_displacement0_vector(jjs*3-1,1);     % 获取当前时刻第jjs个粒子的y方向位移
                    ww_jjs=global_displacement0_vector(jjs*3,1);     % 获取当前时刻第jjs个粒子的z方向位移
                    length1_i_jjs=sqrt((xx_jjs-xx_i)*(xx_jjs-xx_i)+(yy_jjs-yy_i)*(yy_jjs-yy_i)+(zz_jjs-zz_i)*(zz_jjs-zz_i)); % 初始构型下粒子i与粒子jjs之间的距离为length1_i_jjs
                    length2_i_jjs=sqrt((xx_jjs-xx_i+uu_jjs-uu_i)*(xx_jjs-xx_i+uu_jjs-uu_i)+(yy_jjs-yy_i+vv_jjs-vv_i)*(yy_jjs-yy_i+vv_jjs-vv_i)+(zz_jjs-zz_i+ww_jjs-ww_i)*(zz_jjs-zz_i+ww_jjs-ww_i)); % 当前构型下粒子i与粒子jjs之间的距离为length2_i_jjs
                    s1=(length2_i_jjs-length1_i_jjs)/length1_i_jjs;   % 粒子i与粒子jj之间键的伸长率为s0
                    length1_j_iin=sqrt((xx_iin-xx_j)*(xx_iin-xx_j)+(yy_iin-yy_j)*(yy_iin-yy_j)+(zz_iin-zz_j)*(zz_iin-zz_j));  % 初始构型下粒子j与粒子iin之间的距离为length1_j_iin 
                    length2_j_iin=sqrt((xx_iin-xx_j+uu_iin-uu_j)*(xx_iin-xx_j+uu_iin-uu_j)+(yy_iin-yy_j+vv_iin-vv_j)*(yy_iin-yy_j+vv_iin-vv_j)+(zz_iin-zz_j+ww_iin-ww_j)*(zz_iin-zz_j+ww_iin-ww_j));  % 当前构型下粒子j与粒子iin之间的距离为length2_j_iin 
                    s2=(length2_j_iin-length1_j_iin)/length1_j_iin;   % 粒子i与粒子jj之间键的伸长率为s0
                    force_density_x=force_density_x+8*bs*detla_s*(s1-s2)*(xx_jjs-xx_i)/length1_i_jjs; % 粒子i与粒子jjs在x方向的作用加到粒子i所受到的x方向力密度矢量中
                    force_density_y=force_density_y+8*bs*detla_s*(s1-s2)*(yy_jjs-yy_i)/length1_i_jjs; % 粒子i与粒子jjs在y方向的作用加到粒子i所受到的y方向力密度矢量中
                    force_density_z=force_density_z+8*bs*detla_s*(s1-s2)*(zz_jjs-zz_i)/length1_i_jjs; % 粒子i与粒子jjs在z方向的作用加到粒子i所受到的z方向力密度矢量中
                end
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        global_accelerate_vector(i*3-2,1)=(force_density_x*volume_particle+global_load_density_vector(i*3-2,1))/density; % 当前时刻粒子i在x方向的加速度
        global_accelerate_vector(i*3-1,1)=(force_density_y*volume_particle+global_load_density_vector(i*3-1,1))/density; % 当前时刻粒子i在y方向的加速度
        global_accelerate_vector(i*3,1)=(force_density_z*volume_particle+global_load_density_vector(i*3,1))/density; % 当前时刻粒子i在z方向的加速度
end
global_velocity_vector=global_velocity0_vector+global_accelerate_vector*dt;          % 更新当前时刻所有粒子的速度
global_displacement_vector=global_displacement0_vector+global_velocity_vector*dt;          % 更新当前时刻所有粒子的位移

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(itime,output_interval)==0  % 如果itime对output_interval取余数为0，则执行条件语句中的内容
    str1='displacement_';           % 定义输出文件名称中的第1部分 displacement_
    str2=num2str(itime);           % 定义输出文件名称中的第2部分 itime
    str3='.txt';                    % 定义输出文件名称中的第3部分 .txt，即文件拓展名
    str4=strcat(str1,str2,str3);    % 将文件名称写到一个变量中，即displacement_itime.txt，其中itime指具体计算步数，如果itime=1,则输出文件名为displacement_1.txt
    dlmwrite(str4,global_displacement_vector);  % 调用 dlmwrite 函数将global_displacement_vector写入文本文件displacement_itime.txt
end