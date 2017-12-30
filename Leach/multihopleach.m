clear;close all;clc;
con=0;
xm = 200;ym = 200;
sink.x =0.5 * xm;
sink.y = ym + 50;
n =100;    
p=0.05;%(成为簇头的最佳比例)
%Energy Model 
Eo = 1;%(最初能量)
%Eelec=Etx=Erx
ETX1=0.001;
ERX1=0.001;
ETX2=0.0001;
ERX2=0.0001;
%Transmit Amplifier types  Efs,Emp
%Data Aggregation Energy
EDA=5*0.000000001;
INFINITY = 999999999999999;   
%maximum number of rounds
rmax=20;                   
  %循环次数
dead=0;
 
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
 
%随机分布节点
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;%坐标
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    %initially there are no cluster heads only nodes
    S(i).type=0;%普通节点
    S(i).E=Eo;
    plot(S(i).xd,S(i).yd,'bo');  
    hold on;  
end
%  pause;
%计算节点间的距离
for i=1:1:n
    for j=1:1:n
        distance(i,j)=((S(i).xd-(S(j).xd))^2 + (S(i).yd-(S(j).yd) )^2 )^(1\2);
    end
end
 
%覆盖半径
R=((xm*ym)/(5*pi))^(1\2);
%增加通信距离
L=75;
 
%基站坐标
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'x');  
 
%大循环轮数
for r=1:1:45    
     ss=rem((r+3),4)+1; 
     if rem(r-1,4)==0
     figure;        
     end
         E=0;        
    for i=1:1:n              
        E=E+S(i).E;       
        subplot(2,2,ss),plot (S(i).xd,S(i).yd,'bo');
        hold on;
        subplot(2,2,ss),plot (S(n+1).xd,S(n+1).yd,'x');
        if S(i).type==1           
            S(i).type=0;
        end     
    end
%计数死亡节点
  
    for i=1:1:n
        %checking if there is a dead node
        if (S(i).E<=0)
            subplot(2,2,ss),plot(S(i).xd,S(i).yd,'b*');  
            hold on;
            if S(i).type==0%type=0为死亡节点
                dead=dead+1;
                 S(i).type=2;
            end
        end
    end 
	if (n-dead<=10)%节点全部死亡退出循环
        break;
    end
   
    countCHs=1;
    cluster=1; 
       
    for k=1:1:10 %事先决定簇头个数
        for i=1:1:n
            c(i)=rand;
            if(S(i).E*n>=E && S(i).type ==0)
                %如果该节点在候选集合中
                %Election of Cluster Heads
                if( c(i) <= (p/(1-p*mod(r,round(1/p)))))
                    countCHs = countCHs+1;
                    S(i).type = 1;
                    C(cluster).type=S(i).type;
                    C(cluster).xd = S(i).xd;
                    C(cluster).yd = S(i).yd;
                    C(cluster).id = i;
                    C(cluster).E=S(i).E;            
                    dis(cluster) = sqrt( (S(n+1).xd-S(i).xd)^2 + (S(n+1).yd-S(i).yd)^2 );
                    cluster=cluster+1;
                    subplot(2,2,ss), plot(S(i).xd,S(i).yd,'r*');   %簇头节点   
                    X=[S(i).xd,S(n+1).xd]; Y=[S(i).yd,S(n+1).yd];
                    for j=1:1:n
                        if distance(i,j)<=R && j~=i && S(j).type==0
                            S(j).type=3;
                        end
                    end 
                    if countCHs==11
                        con=1;
                        break;
                    end
                end          
            end    
        end  
        if con==1
            con=0;
            break;
        end
    end

    
    
    
    
    
    for i=1:1:n
        if S(i).type ==3
              S(i).type =0;
        end
        min_dis=40000;
        if S(i).type ==0
            for c=1:1:cluster -1 %簇头数量一共是cluster - 1
                temp =sqrt( (S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2 ) ;
                if ( temp < min_dis )
                    min_dis = temp;
                    j=c;
                end     
            end     
            if min_dis<L
                S(i).E = S(i).E - ETX1 * min_dis- ERX1 * min_dis-ETX2 * min_dis;   % 接到簇头加入信息
                C(j).E = C(j).E - ETX1 * min_dis- ERX1 * min_dis-ERX2 * min_dis;   % 簇头发送信息并接收
                S(C(j).id).E=C(j).E; 
                X=[S(i).xd,C(j).xd]; Y=[S(i).yd,C(j).yd];
                subplot(2,2,ss),plot(X,Y);   
                hold on;
            end      
         end
    end
  
	min_dis=40000;
	for i=1:1:cluster -1     
         if dis(i)<min_dis
             min_dis=dis(i);
             c=i;
         end
	end
	X=[S(n+1).xd,C(c).xd]; Y=[S(n+1).yd,C(c).yd];
	subplot(2,2,ss), plot(X,Y,'r-');
	for i=1:1:cluster -1
        if i~=c && abs(C(i).yd-C(c).yd)<=50 && dis(i)<distance(C(i).id,C(c).id)
            X=[S(n+1).xd,C(i).xd]; Y=[S(n+1).yd,C(i).yd]; 
            subplot(2,2,ss), plot(X,Y,'r-');
        end             
	end
  
    for i=1:1:cluster -1
        for j=1:1:cluster -1
            if abs(C(i).xd-C(j).xd)<=70 && abs(C(i).yd-C(j).yd)<=100 %&& distance(C(i).id,C(j).id)<2L
                X=[C(i).xd,C(j).xd]; Y=[C(i).yd,C(j).yd]; 
               subplot(2,2,ss), plot(X,Y,'g-');
               hold on;
            end
        end
    end
        
    for i=1:1:n
        if S(i).type ==1
            S(i).E=S(i).E - ETX2 *( sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 ));
        end
	end
end