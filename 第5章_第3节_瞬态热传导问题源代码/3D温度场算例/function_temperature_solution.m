function[global_temperature_vector]=function_temperature_solution(dimension,size_particle,size_temperature_boundary_left,particle_boundary_left, ...
         size_temperature_boundary_right,particle_boundary_right,size_temperature_boundary_bottom,particle_boundary_bottom,size_temperature_boundary_top, ...
         particle_boundary_top,size_temperature_boundary_behind,particle_boundary_behind,size_temperature_boundary_front,particle_boundary_front, ...
         size_heat_flux_left,particle_heat_flux_left,size_heat_flux_right,particle_heat_flux_right,size_heat_flux_bottom,particle_heat_flux_bottom, ...
         size_heat_flux_top,particle_heat_flux_top,size_heat_flux_behind,particle_heat_flux_behind,size_heat_flux_front,particle_heat_flux_front, ...
         temperature_left,temperature_right,temperature_bottom,temperature_top,temperature_behind,temperature_front,heat_flux_left,heat_flux_right, ...
         heat_flux_bottom,heat_flux_top,heat_flux_behind,heat_flux_front,body_heat_source,boundary_condition_left,boundary_condition_right,boundary_condition_bottom, ...
         boundary_condition_top,boundary_condition_behind,boundary_condition_front,dettx,global_temperature_vector,itime,density,c0,cv,particle,horizon_particle,size_horizon_particle,dt,output_interval)

global_heat_flux_vector=zeros(size_particle,1); % 总体热流矢量赋初值
global_temperature0_vector=global_temperature_vector;  % 将上一时间步计算得到的总体温度矢量赋值给 global_temperature0_vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 一维情况
if dimension==1
%%%%%%%%%%%%%%%%%%%%%%%
% 施加温度边界条件
    % 左边界
    if boundary_condition_left==1                                % 如果左边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_left                     % 对左边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_left(i,1);                           % ii表示左边界施加温度边界条件的粒子编号
            jj=particle_boundary_left(i,2);                           % jj表示与左边界施加温度边界条件的粒子关于左边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_left*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
    % 右边界
    if boundary_condition_right==1                           % 如果右边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_right                % 对右边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_right(i,1);                      % ii表示右边界施加温度边界条件的粒子编号
            jj=particle_boundary_right(i,2);                           % jj表示与右边界施加温度边界条件的粒子关于左边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_right*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加热流边界条件
    % 左边界
    if boundary_condition_left==2    % 如果左边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_left  % 对左边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_left(i);   % ii表示左边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_left/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 右边界
    if boundary_condition_right==2    % 如果右边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_right  % 对右边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_right(i);   % ii表示右边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_right/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加内热源
    if body_heat_source~=0    % 如果施加了内热源，则执行条件语句中的内容
        for i=1:size_particle % 对所有粒子循环
            global_heat_flux_vector(i,1)=global_heat_flux_vector(i,1)+body_heat_source; % 在总体热流矢量对应位置增加热流
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 温度求解
    thermal_response_vector=zeros(size_particle,1);    % 热响应矢量清零
    for i=1:size_particle    % 对所有粒子循环
        xx_i=particle(i,1);     % 获取第i个粒子的x方向坐标
        temperature_i=global_temperature0_vector(i,1);     % 获取当前时刻第i个粒子的温度
        for j=1:size_horizon_particle       % 对第i个粒子作用范围内的所有粒子进行循环
            jj=horizon_particle(i,j);       % 第i个粒子作用范围内的第j个粒子为粒子jj
            if jj>0       % 如果第i个粒子作用范围内的第j个粒子存在，则执行循环条件内的语句
            	xx_j=particle(jj,1);     % 获取第jj个粒子的x方向坐标
                temperature_j=global_temperature0_vector(jj,1);     % 获取当前时刻第jj个粒子的温度
                length1_i_jj=abs(xx_j-xx_i);   % 初始构型下粒子i与粒子jj之间的距离为length_i_jj
                thermal_response_vector(i,1)=thermal_response_vector(i,1)+c0*(temperature_j-temperature_i)/length1_i_jj; % 将粒子jj对粒子i的作用加到热响应矢量中
            end
        end
    end
    global_temperature_vector=global_temperature0_vector+dt*(thermal_response_vector+global_heat_flux_vector)/(density*cv); % 更新当前时刻所有粒子的温度
%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%二维情况
if dimension==2
%%%%%%%%%%%%%%%%%%%%%%%
% 施加温度边界条件
    % 左边界
    if boundary_condition_left==1                                % 如果左边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_left                     % 对左边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_left(i,1);                           % ii表示左边界施加温度边界条件的粒子编号
            jj=particle_boundary_left(i,2);                           % jj表示与左边界施加温度边界条件的粒子关于左边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_left*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
    % 右边界
    if boundary_condition_right==1                           % 如果右边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_right                % 对右边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_right(i,1);                      % ii表示右边界施加温度边界条件的粒子编号
            jj=particle_boundary_right(i,2);                           % jj表示与右边界施加温度边界条件的粒子关于右边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_right*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
    % 下边界
    if boundary_condition_bottom==1                                % 如果下边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_bottom                     % 对下边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_bottom(i,1);                      % ii表示下边界施加温度边界条件的粒子编号
            jj=particle_boundary_bottom(i,2);                           % jj表示与下边界施加温度边界条件的粒子关于下边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_bottom*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
    % 上边界
    if boundary_condition_top==1                                % 如果上边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_top                     % 对上边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_top(i,1);                      % ii表示上边界施加温度边界条件的粒子编号
            jj=particle_boundary_top(i,2);                           % jj表示与上边界施加温度边界条件的粒子关于上边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_top*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加热流边界条件
    % 左边界
    if boundary_condition_left==2    % 如果左边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_left  % 对左边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_left(i);   % ii表示左边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_left/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 右边界
    if boundary_condition_right==2    % 如果右边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_right  % 对右边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_right(i);   % ii表示右边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_right/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 下边界
    if boundary_condition_bottom==2    % 如果下边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_bottom  % 对下边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_bottom(i);   % ii表示下边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_bottom/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 上边界
    if boundary_condition_top==2    % 如果上边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_top  % 对上边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_top(i);   % ii表示上边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_top/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加内热源
    if body_heat_source~=0    % 如果施加了内热源，则执行条件语句中的内容
        for i=1:size_particle % 对所有粒子循环
            global_heat_flux_vector(i,1)=global_heat_flux_vector(i,1)+body_heat_source; % 在总体热流矢量对应位置增加热流
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 温度求解
    thermal_response_vector=zeros(size_particle,1);    % 热响应矢量清零
    for i=1:size_particle    % 对所有粒子循环
        xx_i=particle(i,1);     % 获取第i个粒子的x方向坐标
        yy_i=particle(i,2);     % 获取第i个粒子的y方向坐标
        temperature_i=global_temperature0_vector(i,1);     % 获取当前时刻第i个粒子的温度
        for j=1:size_horizon_particle       % 对第i个粒子作用范围内的所有粒子进行循环
            jj=horizon_particle(i,j);       % 第i个粒子作用范围内的第j个粒子为粒子jj
            if jj>0       % 如果第i个粒子作用范围内的第j个粒子存在，则执行循环条件内的语句
            	xx_j=particle(jj,1);     % 获取第jj个粒子的x方向坐标
                yy_j=particle(jj,2);     % 获取第jj个粒子的y方向坐标
                temperature_j=global_temperature0_vector(jj,1);     % 获取当前时刻第jj个粒子的温度
                length1_i_jj=sqrt((xx_j-xx_i)*(xx_j-xx_i)+(yy_j-yy_i)*(yy_j-yy_i));   % 初始构型下粒子i与粒子jj之间的距离为length_i_jj
                thermal_response_vector(i,1)=thermal_response_vector(i,1)+c0*(temperature_j-temperature_i)/length1_i_jj;  % 将粒子jj对粒子i的作用加到热响应矢量中
            end
        end
    end
    global_temperature_vector=global_temperature0_vector+dt*(thermal_response_vector+global_heat_flux_vector)/(density*cv); % 更新当前时刻所有粒子的温度
%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%三维情况
if dimension==3
%%%%%%%%%%%%%%%%%%%%%%%
% 施加温度边界条件
    % 左边界
    if boundary_condition_left==1                                % 如果左边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_left                     % 对左边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_left(i,1);                           % ii表示左边界施加温度边界条件的粒子编号
            jj=particle_boundary_left(i,2);                           % jj表示与左边界施加温度边界条件的粒子关于左边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_left*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
    % 右边界
    if boundary_condition_right==1                           % 如果右边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_displacement_boundary_right                % 对右边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_right(i,1);                      % ii表示右边界施加温度边界条件的粒子编号
            jj=particle_boundary_right(i,2);                           % jj表示与右边界施加温度边界条件的粒子关于右边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_right*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
    % 下边界
    if boundary_condition_bottom==1                                % 如果下边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_bottom                     % 对下边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_bottom(i,1);                      % ii表示下边界施加温度边界条件的粒子编号
            jj=particle_boundary_bottom(i,2);                           % jj表示与下边界施加温度边界条件的粒子关于下边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_bottom*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
    % 上边界
    if boundary_condition_top==1                                % 如果上边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_top                     % 对上边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_top(i,1);                      % ii表示上边界施加温度边界条件的粒子编号
            jj=particle_boundary_top(i,2);                           % jj表示与上边界施加温度边界条件的粒子关于上边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_top*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
    % 后边界
    if boundary_condition_behind==1                                % 如果后边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_behind                   % 对后边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_behind(i,1);                      % ii表示后边界施加温度边界条件的粒子编号
            jj=particle_boundary_behind(i,2);                           % jj表示与后边界施加温度边界条件的粒子关于后边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_behind*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
    % 前边界
    if boundary_condition_front==1                                % 如果前边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_front                   % 对前边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_front(i,1);                      % ii表示前边界施加温度边界条件的粒子编号
            jj=particle_boundary_front(i,2);                           % jj表示与前边界施加温度边界条件的粒子关于前边界对称的粒子的编号
            global_temperature0_vector(ii,1)=temperature_front*2-global_temperature0_vector(jj,1);     % 在总体温度矢量的对应位置赋予该时间步的温度
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加热流边界条件
    % 左边界
    if boundary_condition_left==2    % 如果左边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_left  % 对左边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_left(i);   % ii表示左边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_left/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 右边界
    if boundary_condition_right==2    % 如果右边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_right  % 对右边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_right(i);   % ii表示右边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_right/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 下边界
    if boundary_condition_bottom==2    % 如果下边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_bottom  % 对下边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_bottom(i);   % ii表示下边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_bottom/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 上边界
    if boundary_condition_top==2    % 如果上边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_top  % 对上边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_top(i);   % ii表示上边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_top/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 后边界
    if boundary_condition_behind==2    % 如果后边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_behind  % 对后边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_behind(i);   % ii表示后边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_behind/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 前边界
    if boundary_condition_front==2    % 如果前边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_front  % 对前边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_front(i);   % ii表示前边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_front/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加内热源
    if body_heat_source~=0    % 如果施加了内热源，则执行条件语句中的内容
        for i=1:size_particle % 对所有粒子循环
            global_heat_flux_vector(i,1)=global_heat_flux_vector(i,1)+body_heat_source; % 在总体热流矢量对应位置增加热流
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 温度求解
    thermal_response_vector=zeros(size_particle,1);    % 热响应矢量清零
    for i=1:size_particle    % 对所有粒子循环
        xx_i=particle(i,1);     % 获取第i个粒子的x方向坐标
        yy_i=particle(i,2);     % 获取第i个粒子的y方向坐标
        zz_i=particle(i,3);     % 获取第i个粒子的y方向坐标
        temperature_i=global_temperature0_vector(i,1);     % 获取当前时刻第i个粒子的温度
        for j=1:size_horizon_particle       % 对第i个粒子作用范围内的所有粒子进行循环
            jj=horizon_particle(i,j);       % 第i个粒子作用范围内的第j个粒子为粒子jj
            if jj>0       % 如果第i个粒子作用范围内的第j个粒子存在，则执行循环条件内的语句
            	xx_j=particle(jj,1);     % 获取第jj个粒子的x方向坐标
                yy_j=particle(jj,2);     % 获取第jj个粒子的y方向坐标
                zz_j=particle(jj,3);     % 获取第jj个粒子的y方向坐标
                temperature_j=global_temperature0_vector(jj,1);     % 获取当前时刻第jj个粒子的温度
                length1_i_jj=sqrt((xx_j-xx_i)*(xx_j-xx_i)+(yy_j-yy_i)*(yy_j-yy_i)+(zz_j-zz_i)*(zz_j-zz_i));   % 初始构型下粒子i与粒子jj之间的距离为length_i_jj
                thermal_response_vector(i,1)=thermal_response_vector(i,1)+c0*(temperature_j-temperature_i)/length1_i_jj;  % 将粒子jj对粒子i的作用加到热响应矢量中
            end
        end
    end
    global_temperature_vector=global_temperature0_vector+dt*(thermal_response_vector+global_heat_flux_vector)/(density*cv); % 更新当前时刻所有粒子的温度
%%%%%%%%%%%%%%%%%%%%%%%
end
if mod(itime,output_interval)==0  % 如果itime对output_interval取余数为0，则执行条件语句中的内容
    str1='temperature_';           % 定义输出文件名称中的第1部分 temperature_
    str2=num2str(itime);           % 定义输出文件名称中的第2部分 itime
    str3='.txt';                    % 定义输出文件名称中的第3部分 .txt，即文件拓展名
    str4=strcat(str1,str2,str3);    % 将文件名称写到一个变量中，即temperature_itime.txt，其中itime指具体计算步数，如果itime=1,则输出文件名为temperature_1.txt
    dlmwrite(str4,global_temperature_vector);  % 调用 dlmwrite 函数将global_temperature_vector写入文本文件temperature_itime.txt
end
