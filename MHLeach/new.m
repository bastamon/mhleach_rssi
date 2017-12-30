clear;%清除变量
xm=100;%设置区域为100*100
ym=100;
sink.x=0.5*xm;%sink（汇聚）节点坐标
sink.y=0.5*ym;
n=100; %区域内的节点数目
p=0.1;% 节点成为簇头的概率 
Eo=0.5;%节点初始能量
ETX=50*0.000000001;%发射单位报文损耗能量 
ERX=50*0.000000001;%接收单位报文损耗能量
Efs=10*0.000000000001;%自由空间能量
Emp=0.0013*0.000000000001;%衰减空间能量
EDA=5*0.000000001;%多路径衰减能量
rmax=3000;%最大的轮数
do=sqrt(Efs/Emp); %计算do 通信半径

R0=20;%自定义的竞争半径的最大值
c=1;

figure(1);%输出图形
for i=1:1:n %i为矩阵1到n，间距为1
S(i).xd=rand(1,1)*xm;%1行1列矩阵
XR(i)=S(i).xd;%随机生成的X轴
S(i).yd=rand(1,1)*ym;
YR(i)=S(i).yd;%随机生成的Y轴
S(i).G=0;%
S(i).type='N';%节点类型为普通
temp_rnd0=i;%随机数值
S(i).E=Eo;%设置初始能量为E0

plot(S(i).xd,S(i).yd,'o');%输出节点，用o表示
hold on;
end
 
S(n+1).xd=sink.x;%汇聚节点X轴坐标
S(n+1).yd=sink.y;%汇聚节点Y轴坐标
plot(S(n+1).xd,S(n+1).yd,'x'); %输出汇聚节点，用x表示 
 
figure(1);
countCHs=0;
rcountCHs=0;
cluster=1;
countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;%第一个节点死亡的标志变量
 
for r=0:1:rmax%r为矩阵0到最大，间距为1
r;

if(mod(r, round(1/p) )==0)%余数为0
for i=1:1:n%i为矩阵1到n，间距为1
S(i).G=0;%簇头数目
S(i).cl=0;
end
end
 


 
 
hold off;
 
dead=0;%节点死亡数
dead_a=0;%高级节点死亡数
dead_n=0;%普通节点死亡数
 
packets_TO_BS=0;%传输sink节点报文数
packets_TO_CH=0;%传输簇头的报文数
PACKETS_TO_CH(r+1)=0;%每轮传送到簇头的报文数 
PACKETS_TO_BS(r+1)=0;%每轮传送到基站的报文数 
 
figure(4);
 
 
for i=1:1:n %i为矩阵1到n，间距为1
if (S(i).E<=0)%检查是否有节点死亡
plot(S(i).xd,S(i).yd,'red .')%输出节点，用红.表示
dead=dead+1;%节点死亡数+1

hold on; 
end
if S(i).E>0%节点能量大于0
S(i).type='N';%节点类型为普通

plot(S(i).xd,S(i).yd,'o');

hold on;
end
end
plot(S(n+1).xd,S(n+1).yd,'x');
 
 
STATISTICS(r+1).DEAD=dead;%r轮后死亡节点数
DEAD(r+1)=dead;%r轮后死亡节点数
DEAD_N(r+1)=dead_n;%r轮后普通节点死亡数
DEAD_A(r+1)=dead_a;%r轮后高级节点死亡数
 
if (dead==1)%第一个节点死亡
if(flag_first_dead==0)%第一个节点死亡周期
first_dead=r;%第一个节点死亡轮数
flag_first_dead=1;%第一个死亡节点标志
end
end
 
 
 
 
 %计算Dmax和Dmin
for i=1:1:n
 S(i).d=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
end
Dmax= S(1).d;
Dmin=S(1).d;
for i=2:1:n
    if(Dmax<=S(i).d)
        Dmax=S(i).d;
    end
    if(Dmin>=S(i).d)
        Dmin=S(i).d;
    end            
end
 
countCHs=0;%簇头的个数
cluster=1;%簇头的数目
for i=1:1:n%i为矩阵1到n，间距为1
if(S(i).E>0)%节点剩余能量大于0
temp_rand=rand; 
if ( (S(i).G)<=0)%没有簇头

    
if( temp_rand <= ( p / ( 1 - p * mod(r,round(1/p)) )) )%普通节点的簇头选举
 %候选簇头
 S(i).G=1;%1表示候选簇首集合
 S(i).Rc=(1-c*(S(i).d-Dmin)/(Dmax-Dmin))*R0;  %%%候选簇首节点的竞争半径

%   if ( S(i).Rc>do)%距离大于通信半径
% S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*(  S(i).Rc* S(i).Rc* S(i).Rc* S(i).Rc)); %能量消耗
% end
% if ( S(i).Rc<=do)%距离小于通信半径
% S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Efs*4000*(  S(i).Rc *  S(i).Rc)); %能量消耗
% end  
 
end 
end
end 
end
%在候选簇头竞争半径范围内选择剩余能量大的为簇头
for i=1:1:n
 if(S(i).E>0)
   if (S(i).G<=0)
          for j=1:1:n     
             if(S(j).E>0)
               if S(j).G==1
                   dist_G = sqrt((S(i).xd-S(j).xd).^2)+((S(i).yd-S(j).yd).^2); %%候选簇首节点间的距离
                 if dist_G< S(j).Rc
                     if S(i).E<S(j).E
                    
                         
                   
countCHs=countCHs+1;%簇头数+1
packets_TO_BS=packets_TO_BS+1;%传送到基站的计数器+1
PACKETS_TO_BS(r+1)=packets_TO_BS;%每轮传送到基站的计数器=传送到基站的计数器
 
S(j).type='C';%节点类型为簇头
S(j).G=100;
C(cluster).xd=S(j).xd;%簇头X轴坐标 
C(cluster).yd=S(j).yd;%簇头Y轴坐标
plot(S(j).xd,S(j).yd,'k*');%输出节点，用黑*表示
 
distance=sqrt( (S(j).xd-(S(n+1).xd) )^2 + (S(j).yd-(S(n+1).yd) )^2 );%计算距离
C(cluster).distance=distance;%距离
C(cluster).id=j;%簇头的节点编号
X(cluster)=S(j).xd;%X轴坐标
Y(cluster)=S(j).yd;%Y轴坐标
cluster=cluster+1;%簇头总数+1
 
 
distance;
if (distance>do)%距离大于通信半径
S(j).E=S(j).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); %能量消耗
end
if (distance<=do)%距离小于通信半径
S(j).E=S(j).E- ( (ETX+EDA)*(4000) + Efs*4000*( distance * distance )); %能量消耗
end  
                         
                         
                         
                         
                         
                         
                     
                     else
                  
        
countCHs=countCHs+1;%簇头数+1
packets_TO_BS=packets_TO_BS+1;%传送到基站的计数器+1
PACKETS_TO_BS(r+1)=packets_TO_BS;%每轮传送到基站的计数器=传送到基站的计数器
 
S(i).type='C';%节点类型为簇头
S(i).G=100;
C(cluster).xd=S(i).xd;%簇头X轴坐标 
C(cluster).yd=S(i).yd;%簇头Y轴坐标
plot(S(i).xd,S(i).yd,'k*');%输出节点，用黑*表示
 
distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%计算距离
C(cluster).distance=distance;%距离
C(cluster).id=i;%簇头的节点编号
X(cluster)=S(i).xd;%X轴坐标
Y(cluster)=S(i).yd;%Y轴坐标
cluster=cluster+1;%簇头总数+1
 
 
distance;
if (distance>do)%距离大于通信半径
S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); %能量消耗
end
if (distance<=do)%距离小于通信半径
S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Efs*4000*( distance * distance )); %能量消耗
end
                     
                     
                     
                     end
                 else
                     %%%%%%%%%%%12/25
                 countCHs=countCHs+1;%簇头数+1
                packets_TO_BS=packets_TO_BS+1;%传送到基站的计数器+1
                PACKETS_TO_BS(r+1)=packets_TO_BS;%每轮传送到基站的计数器=传送到基站的计数器
 
               S(i).type='C';%节点类型为簇头
               S(i).G=100;
              C(cluster).xd=S(i).xd;%簇头X轴坐标 
             C(cluster).yd=S(i).yd;%簇头Y轴坐标
              plot(S(i).xd,S(i).yd,'k*');%输出节点，用黑*表示
 
             distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%计算距离
             C(cluster).distance=distance;%距离
             C(cluster).id=i;%簇头的节点编号
             X(cluster)=S(i).xd;%X轴坐标
             Y(cluster)=S(i).yd;%Y轴坐标
             cluster=cluster+1;%簇头总数+1
 
 
             distance;
             if (distance>do)%距离大于通信半径
             S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); %能量消耗
             end
             if (distance<=do)%距离小于通信半径
             S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Efs*4000*( distance * distance )); %能量消耗
             end
                     
                     
                     
                     
                 end
               end
             end
          end
   end
 end
end


%  %正式簇头 
% countCHs=countCHs+1;%簇头数+1
% packets_TO_BS=packets_TO_BS+1;%传送到基站的计数器+1
% PACKETS_TO_BS(r+1)=packets_TO_BS;%每轮传送到基站的计数器=传送到基站的计数器
%  
% S(i).type='C';%节点类型为簇头
% S(i).G=100;
% C(cluster).xd=S(i).xd;%簇头X轴坐标 
% C(cluster).yd=S(i).yd;%簇头Y轴坐标
% plot(S(i).xd,S(i).yd,'k*');%输出节点，用黑*表示
%  
% distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%计算距离
% C(cluster).distance=distance;%距离
% C(cluster).id=i;%簇头的节点编号
% X(cluster)=S(i).xd;%X轴坐标
% Y(cluster)=S(i).yd;%Y轴坐标
% cluster=cluster+1;%簇头总数+1
%  
%  
% distance;
% if (distance>do)%距离大于通信半径
% S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); %能量消耗
% end
% if (distance<=do)%距离小于通信半径
% S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Efs*4000*( distance * distance )); %能量消耗
% end



 
 
 
STATISTICS(r+1).CLUSTERHEADS=cluster-1;%r轮后簇头数
CLUSTERHS(r+1)=cluster-1;%r轮后簇头数
 
 
for i=1:1:n
if ( S(i).type=='N' && S(i).E>0 )%选举正常节点的相关簇头
if(cluster-1>=1)%簇头总数大于2个
min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );%两节点间最短距离
min_dis_cluster=1;%距离最小的簇头数
for c=1:1:cluster-1
temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
if ( temp<min_dis )
min_dis=temp;
min_dis_cluster=c;
end
end
min_dis;
if (min_dis>do)
S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
end
if (min_dis<=do)
S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
end
 
if(min_dis>0)%能量消散
S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
PACKETS_TO_CH(r+1)=n-dead-cluster+1; 
end
S(i).min_dis=min_dis;
S(i).min_dis_cluster=min_dis_cluster;
end
end
end
hold on;
countCHs;
rcountCHs=rcountCHs+countCHs;
STATISTICS(r+1).ENERGY=0;
for i=1:1:n%当前节点数
if S(i).E > 0%如果节点i剩余能量大于0
STATISTICS(r+1).ENERGY = STATISTICS(r+1).ENERGY +S(i).E;%r轮后节点剩余能量加上节点i的剩余能量
end
end
 


% [vx,vy]=voronoi(X,Y);
% plot(X,Y,'r*',vx,vy,'b-');
% hold on;
% voronoi(X,Y);
% axis([0 xm 0 ym]);
hold off;
end
for i=2:rmax%当前节点数
mylive(i) = n - STATISTICS(i).DEAD;
myenergy(i) = STATISTICS(i).ENERGY;%剩余能量
end
mylive(1)=100;
myenergy(1)=S(1).E+(n-1)*Eo;
figure(2);%输出图形2
hold on;%保持曲线
plot(mylive,'color','r');%用红色输出存活节点数

xlabel('周期数');
ylabel('存活节点');
title('存活节点图');
figure(3);%输出图形3
hold on;%保持曲线
plot(myenergy,'color','r');%用红色输出剩余能量
xlabel('周期数');
ylabel('剩余能量节点');
title('剩余能量图');