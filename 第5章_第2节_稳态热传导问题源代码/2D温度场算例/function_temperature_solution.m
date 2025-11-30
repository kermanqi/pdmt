function[global_temperature_vector]=function_temperature_solution(dimension,size_particle,global_heat_conduction_matrix,size_temperature_boundary_left,particle_boundary_left, ...
         size_temperature_boundary_right,particle_boundary_right,size_temperature_boundary_bottom,particle_boundary_bottom,size_temperature_boundary_top, ...
         particle_boundary_top,size_temperature_boundary_behind,particle_boundary_behind,size_temperature_boundary_front,particle_boundary_front, ...
         size_heat_flux_left,particle_heat_flux_left,size_heat_flux_right,particle_heat_flux_right,size_heat_flux_bottom,particle_heat_flux_bottom, ...
         size_heat_flux_top,particle_heat_flux_top,size_heat_flux_behind,particle_heat_flux_behind,size_heat_flux_front,particle_heat_flux_front, ...
         temperature_left,temperature_right,temperature_bottom,temperature_top,temperature_behind,temperature_front,heat_flux_left,heat_flux_right, ...
         heat_flux_bottom,heat_flux_top,heat_flux_behind,heat_flux_front,body_heat_source,boundary_condition_left,boundary_condition_right,boundary_condition_bottom, ...
         boundary_condition_top,boundary_condition_behind,boundary_condition_front,volume_particle,dettx)

global_heat_flux_vector=zeros(size_particle,1); % 总体热流矢量赋初值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 一维情况
if dimension==1
%%%%%%%%%%%%%%%%%%%%%%%
% 施加温度边界条件
    % 左边界
    if boundary_condition_left==1                                % 如果左边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_left                     % 对左边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_left(i);                           % ii表示左边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_left;    % 采用乘大数法对总体热流矢量进行操作
        end
    end
    % 右边界
    if boundary_condition_right==1                           % 如果右边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_right                % 对右边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_right(i);                      % ii表示右边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_right; % 采用乘大数法对总体热流矢量进行操作
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加热流边界条件
    % 左边界
    if boundary_condition_left==2    % 如果左边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_left  % 对左边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_left(i);   % ii表示左边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_left*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 右边界
    if boundary_condition_right==2    % 如果右边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_right  % 对右边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_right(i);   % ii表示右边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_right*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加内热源
    if body_heat_source~=0    % 如果施加了内热源，则执行条件语句中的内容
        for i=1:size_particle % 对所有粒子循环
            global_heat_flux_vector(i,1)=global_heat_flux_vector(i,1)+body_heat_source(1)*volume_particle; % 在总体热流矢量对应位置增加热流
        end
    end
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
            ii=particle_boundary_left(i);                           % ii表示左边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_left;    % 采用乘大数法对总体热流矢量进行操作
        end
    end
    % 右边界
    if boundary_condition_right==1                           % 如果右边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_right                % 对右边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_right(i);                      % ii表示右边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_right; % 采用乘大数法对总体热流矢量进行操作
        end
    end
    % 下边界
    if boundary_condition_bottom==1                                % 如果下边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_bottom                     % 对下边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);                           % ii表示下边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_bottom;    % 采用乘大数法对总体热流矢量进行操作
        end
    end
    % 上边界
    if boundary_condition_top==1                                % 如果上边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_top                     % 对上边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_top(i);                           % ii表示上边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_top;    % 采用乘大数法对总体热流矢量进行操作
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加热流边界条件
    % 左边界
    if boundary_condition_left==2    % 如果左边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_left  % 对左边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_left(i);   % ii表示左边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_left*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 右边界
    if boundary_condition_right==2    % 如果右边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_right  % 对右边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_right(i);   % ii表示右边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_right*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 下边界
    if boundary_condition_bottom==2    % 如果下边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_bottom  % 对下边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_bottom(i);   % ii表示下边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_bottom*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 上边界
    if boundary_condition_top==2    % 如果上边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_top  % 对上边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_top(i);   % ii表示上边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_top*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加内热源
    if body_heat_source~=0    % 如果施加了内热源，则执行条件语句中的内容
        for i=1:size_particle % 对所有粒子循环
            global_heat_flux_vector(i,1)=global_heat_flux_vector(i,1)+body_heat_source(1)*volume_particle; % 在总体热流矢量对应位置增加热流
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%三维情况
if dimension==3
%%%%%%%%%%%%%%%%%%%%%%%
% 施加温度边界条件
    % 左边界
    if boundary_condition_left==1                                % 如果左边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_left                     % 对左边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_left(i);                           % ii表示左边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_left;    % 采用乘大数法对总体热流矢量进行操作
        end
    end
    % 右边界
    if boundary_condition_right==1                           % 如果右边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_right                % 对右边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_right(i);                      % ii表示右边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_right; % 采用乘大数法对总体热流矢量进行操作
        end
    end
    % 下边界
    if boundary_condition_bottom==1                                % 如果下边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_bottom                     % 对下边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_bottom(i);                           % ii表示下边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_bottom;    % 采用乘大数法对总体热流矢量进行操作
        end
    end
    % 上边界
    if boundary_condition_top==1                                % 如果上边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_top                     % 对上边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_top(i);                           % ii表示上边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_top;    % 采用乘大数法对总体热流矢量进行操作
        end
    end
    % 后边界
    if boundary_condition_behind==1                                % 如果后边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_behind                   % 对后边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_behind(i);                           % ii表示后边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_behind;    % 采用乘大数法对总体热流矢量进行操作
        end
    end
    % 前边界
    if boundary_condition_front==1                                % 如果前边界定义了温度边界条件，则执行条件语句中的内容
        for i=1:size_temperature_boundary_front                   % 对前边界施加温度边界条件的所有粒子循环
            ii=particle_boundary_front(i);                           % ii表示前边界施加温度边界条件的粒子编号
            global_heat_conduction_matrix(ii,ii)=global_heat_conduction_matrix(ii,ii)*1.0e8;    % 采用乘大数法对总体热传导矩阵进行操作
            global_heat_flux_vector(ii,1)=global_heat_conduction_matrix(ii,ii)*temperature_front;    % 采用乘大数法对总体热流矢量进行操作
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加热流边界条件
    % 左边界
    if boundary_condition_left==2    % 如果左边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_left  % 对左边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_left(i);   % ii表示左边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_left*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 右边界
    if boundary_condition_right==2    % 如果右边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_right  % 对右边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_right(i);   % ii表示右边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_right*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 下边界
    if boundary_condition_bottom==2    % 如果下边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_bottom  % 对下边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_bottom(i);   % ii表示下边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_bottom*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 上边界
    if boundary_condition_top==2    % 如果上边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_top  % 对上边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_top(i);   % ii表示上边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_top*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 后边界
    if boundary_condition_behind==2    % 如果后边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_behind  % 对后边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_behind(i);   % ii表示后边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_behind*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
    % 前边界
    if boundary_condition_front==2    % 如果前边界定义了热流边界条件，则执行条件语句中的内容
        for i=1:size_heat_flux_front  % 对前边界施加热流边界条件的所有粒子循环
            ii=particle_heat_flux_front(i);   % ii表示前边界施加热流边界条件的粒子编号
            global_heat_flux_vector(ii,1)=global_heat_flux_vector(ii,1)+heat_flux_front*volume_particle/dettx; % 在总体热流矢量对应位置增加热流
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
% 施加内热源
    if body_heat_source~=0    % 如果施加了内热源，则执行条件语句中的内容
        for i=1:size_particle % 对所有粒子循环
            global_heat_flux_vector(i,1)=global_heat_flux_vector(i,1)+body_heat_source(1)*volume_particle; % 在总体热流矢量对应位置增加热流
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global_temperature_vector=global_heat_conduction_matrix\global_heat_flux_vector;  % 对总体温度矢量进行求解
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%