function[global_heat_conduction_matrix]=function_global_heat_conduction_matrix(particle,size_particle,horizon_particle,size_horizon_particle,c0_temperature,volume_particle,dimension)


size_vector=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%一维情况
if dimension==1
    for i=1:size_particle                       % 对粒子循环
        xx_i=particle(i,1);                     % 获取第i个粒子的x方向坐标
        for j=1:size_horizon_particle           % 对第i个粒子域内的粒子进行循环
            jj=horizon_particle(i,j);           % 第i个粒子域内的第j个粒子的粒子编号为jj
            if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
                xx_j=particle(jj,1);            % 获取第jj个粒子的x方向坐标
                length_i_jj=abs(xx_j-xx_i);   % 粒子i与粒子j之间的距离为length_i_jj
                kk_PD=c0_temperature/length_i_jj*[1 -1;-1 1]*volume_particle*volume_particle; % 粒子i与粒子j组成键的热传导矩阵
                i1=[i];                 % 需要组装到总体热传导矩阵中的行向量
                j1=[i,jj];     % 需要组装到总体热传导矩阵中的列向量
                for ii3=1:1                     % 对行向量进行循环
                    for jj3=1:2                 % 对列向量进行循环
                        size_vector=size_vector+1;  % 变量+1
                        vector_x(size_vector)=i1(ii3);  % 行向量元素赋值
                        vector_y(size_vector)=j1(jj3);  % 列向量元素赋值
                        heat_conduction_matrix(size_vector)=kk_PD(ii3,jj3);  % 总体热传导矩阵元素赋值
                    end
                end
            end
        end
    end
    global_heat_conduction_matrix=sparse(vector_x,vector_y,heat_conduction_matrix,size_particle,size_particle); % 合材总体热传导矩阵
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%二维情况
if dimension==2
    for i=1:size_particle                       % 对粒子循环
        xx_i=particle(i,1);                     % 获取第i个粒子的x方向坐标
        yy_i=particle(i,2);                     % 获取第i个粒子的y方向坐标
        for j=1:size_horizon_particle           % 对第i个粒子域内的粒子进行循环
            jj=horizon_particle(i,j);           % 第i个粒子域内的第j个粒子的粒子编号为jj
            if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
                xx_j=particle(jj,1);            % 获取第jj个粒子的x方向坐标
                yy_j=particle(jj,2);            % 获取第jj个粒子的y方向坐标
                length_i_jj=sqrt((xx_j-xx_i)*(xx_j-xx_i)+(yy_j-yy_i)*(yy_j-yy_i));   % 粒子i与粒子j之间的距离为length_i_jj
                kk_PD=c0_temperature/length_i_jj*[1 -1;-1 1]*volume_particle*volume_particle; % 粒子i与粒子j组成键的热传导矩阵
                i1=[i];                 % 需要组装到总体热传导矩阵中的行向量
                j1=[i,jj];     % 需要组装到总体热传导矩阵中的列向量
                for ii3=1:1                     % 对行向量进行循环
                    for jj3=1:2                 % 对列向量进行循环
                        size_vector=size_vector+1;  % 变量+1
                        vector_x(size_vector)=i1(ii3);  % 行向量元素赋值
                        vector_y(size_vector)=j1(jj3);  % 列向量元素赋值
                        heat_conduction_matrix(size_vector)=kk_PD(ii3,jj3);  % 总体热传导矩阵元素赋值
                    end
                end
            end
        end
    end
    global_heat_conduction_matrix=sparse(vector_x,vector_y,heat_conduction_matrix,size_particle,size_particle); % 合材总体热传导矩阵
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%三维情况
if dimension==3
    for i=1:size_particle                       % 对粒子循环
        xx_i=particle(i,1);                     % 获取第i个粒子的x方向坐标
        yy_i=particle(i,2);                     % 获取第i个粒子的y方向坐标
        zz_i=particle(i,3);                     % 获取第i个粒子的z方向坐标
        for j=1:size_horizon_particle           % 对第i个粒子域内的粒子进行循环
            jj=horizon_particle(i,j);           % 第i个粒子域内的第j个粒子的粒子编号为jj
            if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
                xx_j=particle(jj,1);            % 获取第jj个粒子的x方向坐标
                yy_j=particle(jj,2);            % 获取第jj个粒子的y方向坐标
                zz_j=particle(jj,3);            % 获取第jj个粒子的z方向坐标
                length_i_jj=sqrt((xx_j-xx_i)*(xx_j-xx_i)+(yy_j-yy_i)*(yy_j-yy_i)+(zz_j-zz_i)*(zz_j-zz_i));   % 粒子i与粒子j之间的距离为length_i_jj
                kk_PD=c0_temperature/length_i_jj*[1 -1;-1 1]*volume_particle*volume_particle; % 粒子i与粒子j组成键的热传导矩阵
                i1=[i];                 % 需要组装到总体热传导矩阵中的行向量
                j1=[i,jj];     % 需要组装到总体热传导矩阵中的列向量
                for ii3=1:1                     % 对行向量进行循环
                    for jj3=1:2                 % 对列向量进行循环
                        size_vector=size_vector+1;  % 变量+1
                        vector_x(size_vector)=i1(ii3);  % 行向量元素赋值
                        vector_y(size_vector)=j1(jj3);  % 列向量元素赋值
                        heat_conduction_matrix(size_vector)=kk_PD(ii3,jj3);  % 总体热传导矩阵元素赋值
                    end
                end
            end
        end
    end
    global_heat_conduction_matrix=sparse(vector_x,vector_y,heat_conduction_matrix,size_particle,size_particle); % 合材总体热传导矩阵
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%