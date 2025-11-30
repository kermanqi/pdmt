clear
rand('state',sum(100*clock));
% 0.35, 0.45, 0.56
tijifenshu=0.3;
radius=0.202;
dettx=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%节点坐标
size_node=256*250;
ii=0;
node=zeros(size_node,2);
for i=1:256
    for j=1:250
        ii=ii+1;
        node(ii,1)=(i-128-0.5)*dettx;
        node(ii,2)=(j-125-0.5)*dettx;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%找圆心
i1=1;
for i=1:1000000
    coords_circle=(rand(1,2)-0.5)*3;
    xx_circle=coords_circle(1,1);
    yy_circle=coords_circle(1,2);
    if i==1
        set_circle(i,1)=xx_circle;
        set_circle(i,2)=yy_circle;
    end
    if i>1
        jj=0;
        for ii=1:(i1)
            xx1_circle=set_circle(ii,1);
            yy1_circle=set_circle(ii,2);
            length=sqrt((xx1_circle-xx_circle)*(xx1_circle-xx_circle)+(yy1_circle-yy_circle)*(yy1_circle-yy_circle));
            if length<radius*2.2
                jj=1;
            end
        end
        if jj<0.5
            i1=i1+1;
            length;
            set_circle(i1,1)=xx_circle;
            set_circle(i1,2)=yy_circle;
        end
    end
    
%节点属性
class_node=zeros(size_node,1);
for i=1:size_node
    xxi=node(i,1);
    yyi=node(i,2);
    for j=1:i1
        xxj=set_circle(j,1);
        yyj=set_circle(j,2);
        length=sqrt((xxj-xxi)*(xxj-xxi)+(yyj-yyi)*(yyj-yyi));
        if length<radius
            class_node(i,1)=1;
            break;
        end
    end
end
tiji=sum(class_node)/size_node
    if tiji>(tijifenshu-0.005)
        break;
    end
end
size_circle=i1;
% scatter(set_circle(:,1),set_circle(:,2))

class2_node=zeros(250,256);
ii=0;
for j=1:256
    for i=1:250
        ii=ii+1;
        class2_node(i,j)=class_node(ii,1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pcolor(class2_node);shading interp;
    
    
    
    
    
    
    
    
    
    
    
    
    
