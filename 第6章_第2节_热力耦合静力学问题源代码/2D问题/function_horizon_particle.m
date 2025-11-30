function[horizon_particle]=function_horizon_particle(particle,size_particle,size_horizon_particle,detla,dimension)

horizon_particle=zeros(size_particle,size_horizon_particle);    %horizon_particle矩阵清零
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%一维情况
if dimension==1    
    for i=1:size_particle %对粒子循环
        xx1=particle(i,1); %获取第i个粒子的坐标
        i1=0; %对第i个粒子域内包含的粒子从0开始计数
        for j=1:size_particle %对粒子循环
            if i~=j %确保第i个粒子域内包含的粒子不是该粒子本身
                xx2=particle(j,1); %获取第j个粒子的坐标
                l=abs(xx1-xx2); %计算第i个粒子与第j个粒子的距离
                if l<=detla	%判断如果两粒子之间距离小于等于粒子作用范围，则第i个粒子域内包含粒子j
                    i1=i1+1; %第i个粒子域内包含的粒子数加1
                    horizon_particle(i,i1)=j; %第i个粒子域内包含的第i1个粒子为粒子j
                end
            end
            if i1==size_horizon_particle
                break; %如果i1等于粒子域内包含的最大粒子个数，则第i个粒子域内包含的粒子已经全部找到，此时结束对第i个粒子的循环
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%二维情况
if dimension==2
    for i=1:size_particle %对粒子循环
        %获取第i个粒子的坐标
        xx1=particle(i,1);
        yy1=particle(i,2);
        i1=0; %对第i个粒子域内包含的粒子从0开始计数
        for j=1:size_particle %对粒子循环
            if i~=j %确保第i个粒子域内包含的粒子不是该粒子本身
                %获取第j个粒子的坐标
                xx2=particle(j,1);
                yy2=particle(j,2);
                l=sqrt((xx1-xx2)^2+(yy1-yy2)^2); %计算第i个粒子与第j个粒子的距离
                if l<=detla %判断如果两粒子之间距离小于等于粒子作用范围，则第i个粒子域内包含粒子j
                    i1=i1+1; %第i个粒子域内包含的粒子数加1
                    horizon_particle(i,i1)=j; %第i个粒子域内包含的第i1个粒子为粒子j
                end
            end
            if i1==size_horizon_particle
                break; %如果i1等于粒子域内包含的最大粒子个数，则第i个粒子域内包含的粒子已经全部找到，此时结束对第i个粒子的循环
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%三维情况
if dimension==3
    for i=1:size_particle %对粒子循环
        %获取第i个粒子的坐标
        xx1=particle(i,1);
        yy1=particle(i,2);
        zz1=particle(i,3);
        i1=0; %对第i个粒子域内包含的粒子从0开始计数
        for j=1:size_particle %对粒子循环
            if i~=j %确保第i个粒子域内包含的粒子不是该粒子本身
                %获取第j个粒子的坐标
                xx2=particle(j,1);
                yy2=particle(j,2);
                zz2=particle(j,3);
                l=sqrt((xx1-xx2)^2+(yy1-yy2)^2+(zz1-zz2)^2); %计算第i个粒子与第j个粒子的距离
                if l<=detla %判断如果两粒子之间距离小于等于粒子作用范围，则第i个粒子域内包含粒子j
                    i1=i1+1; %第i个粒子域内包含的粒子数加1
                    horizon_particle(i,i1)=j; %第i个粒子域内包含的第i1个粒子为粒子j
                end
            end
            if i1==size_horizon_particle %如果i1等于粒子域内包含的最大粒子个数，则第i个粒子域内包含的粒子已经全部找到，此时结束对第i个粒子的循环
                break;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%