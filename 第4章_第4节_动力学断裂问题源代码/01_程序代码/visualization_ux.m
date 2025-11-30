size_x=106;
size_y=100;
global_displacement_vector=load('displacement_9800.txt');
displacement_x=zeros(size_y,size_x);    % x方向位移写为矩阵形式
ii=0;                                   % 变量初始值取0
for i=1:size_y                          % 按y方向粒子总数循环
	for j=1:size_x                      % 按x方向粒子总数循环
        ii=ii+1;                        % ii取值+1，从而遍历每个粒子
        displacement_x(i,j)=global_displacement_vector(ii*2-1,1);   % 总体位移向量中的x方向位移写入x方向位移矩阵
	end
end
pcolor(displacement_x);shading interp;  % x方向位移结果可视化
set(colorbar,'FontSize',14,'FontName','Times New Roman');  % 设置图例的字号字体
colormap Jet  % 图像颜色映射设置为Jet