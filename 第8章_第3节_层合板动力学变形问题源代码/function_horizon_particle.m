function[horizon_particle,horizon_particle_n,horizon_particle_s_lower,horizon_particle_s_upper]= ...
function_horizon_particle(particle,size_particle,size_horizon_particle,detla,size_horizon_particle_n,size_horizon_particle_s,dettx,thick)

horizon_particle=zeros(size_particle,size_horizon_particle);    % 面内作用粒子的矩阵清零
horizon_particle_n=zeros(size_particle,size_horizon_particle_n);    % 面外法向作用粒子的矩阵清零
horizon_particle_s_lower=zeros(size_particle,size_horizon_particle_s);    % 层间剪切作用且位于粒子下层粒子的矩阵清零
horizon_particle_s_upper=zeros(size_particle,size_horizon_particle_s);    % 层间剪切作用且位于粒子上层粒子的矩阵清零
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 寻找与粒子进行面内作用的粒子
for i=1:size_particle % 对所有粒子进行循环
	xx_i=particle(i,1);  % 获取第i个粒子的x坐标
	yy_i=particle(i,2);  % 获取第i个粒子的y坐标
	zz_i=particle(i,3);  % 获取第i个粒子的z坐标
	i1=0; % 对第i个粒子域内包含的粒子从0开始计数
	for j=1:size_particle % 对所有粒子进行循环
        if i~=j  % 如果第i个粒子域内包含的粒子不是该粒子本身，则执行条件语句中的内容
            xx_j=particle(j,1);  % 获取第j个粒子的x坐标
            yy_j=particle(j,2);  % 获取第j个粒子的y坐标
            zz_j=particle(j,3);  % 获取第j个粒子的z坐标
            l=sqrt((xx_i-xx_j)^2+(yy_i-yy_j)^2); % 计算第i个粒子与第j个粒子的距离
            if l<=detla % 如果两粒子之间距离小于等于粒子作用范围，则执行条件语句中的内容
                if abs(zz_j-zz_i)<0.5*thick % 如果两粒子z坐标差的绝对值小于0.5*thick，则两粒子位于同一层，执行条件语句中的内容
                    i1=i1+1; % 第i个粒子域内包含的粒子数加1
                    horizon_particle(i,i1)=j; % 第i个粒子域内包含的第i1个粒子为粒子j
                end
            end
        end
        if i1==size_horizon_particle
            break; % 如果i1等于粒子域内包含的最大粒子个数，则与粒子i进行面内作用的粒子已经全部找到，此时结束对第i个粒子的循环
        end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 寻找与粒子进行面外法向作用的粒子
for i=1:size_particle % 对所有粒子进行循环
	xx_i=particle(i,1);  % 获取第i个粒子的x坐标
	yy_i=particle(i,2);  % 获取第i个粒子的y坐标
	zz_i=particle(i,3);  % 获取第i个粒子的z坐标
	for j=1:size_particle % 对所有粒子进行循环
        if i~=j % 确保第i个粒子域内包含的粒子不是该粒子本身
            xx_j=particle(j,1);  % 获取第j个粒子的x坐标
            yy_j=particle(j,2);  % 获取第j个粒子的y坐标
            zz_j=particle(j,3);  % 获取第j个粒子的z坐标
            if abs(xx_j-xx_i)<0.5*dettx && abs(yy_j-yy_i)<0.5*dettx % 如果两粒子x坐标差的绝对值和y坐标差的绝对值均小于0.5*dettx，则两粒子的x坐标相等，y坐标相等，执行条件语句中的内容
                if (zz_i-zz_j)>0.5*thick && (zz_i-zz_j)<1.5*thick % 如果粒子i的z坐标与粒子j的z坐标的差大于0.5*thick且小于1.5*thick，则粒子j位于粒子i下面一层，执行条件语句中的内容
                    horizon_particle_n(i,1)=j;  % 将粒子j的编号赋予horizon_particle_n(i,1)
                end
                if (zz_j-zz_i)>0.5*thick && (zz_j-zz_i)<1.5*thick % 如果粒子j的z坐标与粒子i的z坐标的差大于0.5*thick且小于1.5*thick，则粒子j位于粒子i上面一层，执行条件语句中的内容
                    horizon_particle_n(i,2)=j; % 将粒子j的编号赋予horizon_particle_n(i,2)
                end
            end
        end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 寻找与粒子进行层间剪切作用的粒子
for i=1:size_particle  % 对所有粒子进行循环
	zz_i=particle(i,3);  % 获取第i个粒子的z坐标
	for j=1:size_horizon_particle  % 对粒子i作用范围内的所有粒子进行循环
        jj=horizon_particle(i,j);           % 第i个粒子域内的第j个粒子的粒子编号为jj
        if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
            xx_jj=particle(jj,1);            % 获取第jj个粒子的x方向坐标
            yy_jj=particle(jj,2);            % 获取第jj个粒子的y方向坐标
            ii_lower=0; % 与粒子i存在层间剪切作用且位于粒子下层的粒子数清零
            ii_upper=0; % 与粒子i存在层间剪切作用且位于粒子上层的粒子数清零
            for k=1:size_particle  % 对所有粒子进行循环
                xx_k=particle(k,1);	% 获取第k个粒子的x坐标
                yy_k=particle(k,2);	% 获取第k个粒子的y坐标
                zz_k=particle(k,3);	% 获取第k个粒子的z坐标
                if abs(xx_k-xx_jj)<0.5*dettx && abs(yy_k-yy_jj)<0.5*dettx % 如果粒子jj和粒子k的x坐标差的绝对值和y坐标差的绝对值均小于0.5*dettx，则两粒子的x坐标相等，y坐标相等，执行条件语句中的内容
                    if (zz_i-zz_k)>0.5*thick && (zz_i-zz_k)<1.5*thick % 如果粒子i的z坐标与粒子k的z坐标的差大于0.5*thick且小于1.5*thick，则粒子k位于粒子i下面一层，执行条件语句中的内容
                        horizon_particle_s_lower(i,j)=k; % 将粒子k的编号赋予horizon_particle_s_lower(i,j)
                        ii_lower=1; % ii_lower赋值为1，表示已经找到与粒子jj的x坐标相等，y坐标相等，且位于粒子i下面一层与粒子i进行层间剪切作用的粒子
                    end
                    if (zz_k-zz_i)>0.5*thick && (zz_k-zz_i)<1.5*thick % 如果粒子k的z坐标与粒子i的z坐标的差大于0.5*thick且小于1.5*thick，则粒子k位于粒子i上面一层，执行条件语句中的内容
                        horizon_particle_s_upper(i,j)=k; % 将粒子k的编号赋予horizon_particle_s_upper(i,j)
                        ii_upper=1; % ii_lower赋值为1，表示已经找到与粒子jj的x坐标相等，y坐标相等，且位于粒子i上面一层与粒子i进行层间剪切作用的粒子
                    end
                    if ii_lower==1 && ii_upper==1
                        break; % 如果ii_lower和ii_upper均为1，表示已经找到与粒子jj的x坐标相等，y坐标相等，且与粒子i进行层间剪切作用的所有粒子，此时结束对粒子i的循环
                    end
                end
            end
        end
	end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%