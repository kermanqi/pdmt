function[classify_bond]=function_classify_bond(particle,size_particle,size_horizon_particle,horizon_particle,theta,dettx)

classify_bond=zeros(size_particle,size_horizon_particle);    % classify_bond矩阵清零
for i=1:size_particle	% 对粒子循环
	xx_i=particle(i,1);	% 获取第i个粒子的x方向坐标
	yy_i=particle(i,2);	% 获取第i个粒子的y方向坐标
	for j=1:size_horizon_particle	% 对第i个粒子域内的粒子进行循环
        jj=horizon_particle(i,j);   % 第i个粒子域内的第j个粒子的粒子编号为jj
        if jj>0 % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
            xx_j=particle(jj,1); % 获取第jj个粒子的x方向坐标
            yy_j=particle(jj,2); % 获取第jj个粒子的y方向坐标
            if theta==0 % 如果单层板纤维与x轴之间夹角为0，则执行条件语句中的内容
                if abs(xx_j-xx_i)>0.1*dettx && abs(yy_j-yy_i)<0.1*dettx % 如果(xx_j-xx_i)的绝对值大于0.1*dettx且(yy_j-yy_i)的绝对值小于0.1*dettx，则执行条件语句中的内容
                    classify_bond(i,j)=1; % 第i个粒子域内的第j个键包含纤维键
                end
            end
            if theta==pi/4 % 如果单层板纤维与x轴之间夹角为pi/4，则执行条件语句中的内容
                if abs((xx_j-xx_i)-(yy_j-yy_i))<0.1*dettx % 如果((xx_j-xx_i)-(yy_j-yy_i))的绝对值，则执行条件语句中的内容
                    classify_bond(i,j)=1; % 第i个粒子域内的第j个键包含纤维键
                end
            end
            if theta==pi/2 % 如果单层板纤维与x轴之间夹角为pi/2，则执行条件语句中的内容
                if abs(xx_j-xx_i)<0.1*dettx && abs(yy_j-yy_i)>0.1*dettx % 如果(xx_j-xx_i)的绝对值小于0.1*dettx且(yy_j-yy_i)的绝对值大于0.1*dettx，则执行条件语句中的内容
                    classify_bond(i,j)=1; % 第i个粒子域内的第j个键包含纤维键
                end
            end
        end
	end
end