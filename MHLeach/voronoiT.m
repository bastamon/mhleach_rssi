
clear;%�������
xm=100;%��������Ϊ100*100
ym=100;
sink.x=0.5*xm;%sink����ۣ��ڵ�����
sink.y=0.5*ym;
n=100; %�����ڵĽڵ���Ŀ
p=0.1;% �ڵ��Ϊ��ͷ�ĸ��� 
Eo=0.5;%�ڵ��ʼ����
ETX=50*0.000000001;%���䵥λ����������� 
ERX=50*0.000000001;%���յ�λ�����������
Efs=10*0.000000000001;%���ɿռ�����
Emp=0.0013*0.000000000001;%˥���ռ�����
EDA=5*0.000000001;%��·��˥������
m=0.1;%��Ϊ�߼��ڵ����
a=1;%����
rmax=1500;%��������
do=sqrt(Efs/Emp); %����do ͨ�Ű뾶
S=struct('xd',zeros(1,1),'yd',zeros(1,1),'G',zeros(1,1),'type',zeros(1,1),'E',zeros(1,1));
XR=zeros(1,100);
YR=zeros(1,100);

figure(1);%���ͼ��
for i=1:1:n %iΪ����1��n�����Ϊ1
    S(i).xd=rand(1,1)*xm;%1��1�о���
    XR(i)=S(i).xd;%������ɵ�X��
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;%������ɵ�Y��
    S(i).G=0;%
    S(i).type='N';%�ڵ�����Ϊ��ͨ
    temp_rnd0=i;%�����ֵ
    S(i).E=Eo;%���ó�ʼ����ΪE0

    plot(S(i).xd,S(i).yd,'o');%����ڵ㣬��o��ʾ
    hold on;
end
pause;
S(n+1).xd=sink.x;%��۽ڵ�X������
S(n+1).yd=sink.y;%��۽ڵ�Y������
plot(S(n+1).xd,S(n+1).yd,'x'); %�����۽ڵ㣬��x��ʾ 
 
figure(1);
countCHs=0;
rcountCHs=0;
cluster=1;
countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;%��һ���ڵ������ı�־����
%%%%%%%%%init
PACKETS_TO_CH=zeros(1,rmax+1);
PACKETS_TO_BS=zeros(1,rmax+1);
STATISTICS=zeros(1,rmax+1);
DEAD=zeros(1,rmax+1);
DEAD_N=zeros(1,rmax+1);
DEAD_A=zeros(1,rmax+1);
for r=0:1:rmax%rΪ����0����󣬼��Ϊ1
    r;

    if(mod(r, round(1/p) )==0)%����Ϊ0
        for i=1:1:n%iΪ����1��n�����Ϊ1
        	S(i).G=0;%��ͷ��Ŀ
            S(i).cl=0;
        end
    end
 


 
 
    hold off;
 
    dead=0;%�ڵ�������
	dead_a=0;%�߼��ڵ�������
	dead_n=0;%��ͨ�ڵ�������
	 
	packets_TO_BS=0;%����sink�ڵ㱨����
	packets_TO_CH=0;%�����ͷ�ı�����
	PACKETS_TO_CH(r+1)=0;%ÿ�ִ��͵���ͷ�ı����� 
	PACKETS_TO_BS(r+1)=0;%ÿ�ִ��͵���վ�ı����� 
	 
	figure(4);
	 
	 
	for i=1:1:n %iΪ����1��n�����Ϊ1
		if (S(i).E<=0)%����Ƿ��нڵ�����
			plot(S(i).xd,S(i).yd,'red .')%����ڵ㣬�ú�.��ʾ
			dead=dead+1;%�ڵ�������+1

			hold on; 
		end
		if S(i).E>0%�ڵ���������0
			S(i).type='N';%�ڵ�����Ϊ��ͨ

			plot(S(i).xd,S(i).yd,'o');

			hold on;
		end
	end
	plot(S(n+1).xd,S(n+1).yd,'x');
	 
	 
	STATISTICS(r+1).DEAD=dead;%r�ֺ������ڵ���
	DEAD(r+1)=dead;%r�ֺ������ڵ���
	DEAD_N(r+1)=dead_n;%r�ֺ���ͨ�ڵ�������
	DEAD_A(r+1)=dead_a;%r�ֺ�߼��ڵ�������
	 
	if (dead==1)%��һ���ڵ�����
		if(flag_first_dead==0)%��һ���ڵ���������
			first_dead=r;%��һ���ڵ���������
			flag_first_dead=1;%��һ�������ڵ��־
		end
	end
	 
	 
	 
	 
	 
	 
	countCHs=0;%��ͷ�ĸ���
	cluster=1;%��ͷ����Ŀ
    C=struct('xd',{},'yd',{});
    
    
	for i=1:1:n%iΪ����1��n�����Ϊ1
		if(S(i).E>0)%�ڵ�ʣ����������0
			temp_rand=rand; 
			if ( (S(i).G)<=0)%û�д�ͷ
			 
				if( temp_rand <= ( p / ( 1 - p * mod(r,round(1/p)) )) )%��ͨ�ڵ�Ĵ�ͷѡ��
				 
					countCHs=countCHs+1;%��ͷ��+1
					packets_TO_BS=packets_TO_BS+1;%���͵���վ�ļ�����+1
					PACKETS_TO_BS(r+1)=packets_TO_BS;%ÿ�ִ��͵���վ�ļ�����=���͵���վ�ļ�����
					 
					S(i).type='C';%�ڵ�����Ϊ��ͷ
					S(i).G=100;
					C(cluster).xd=S(i).xd;%��ͷX������ 
					C(cluster).yd=S(i).yd;%��ͷY������
					plot(S(i).xd,S(i).yd,'k*');%����ڵ㣬�ú�*��ʾ
					 
					distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%�������
					C(cluster).distance=distance;%����
					C(cluster).id=i;%��ͷ�Ľڵ���
					X(cluster)=S(i).xd;%X������
					Y(cluster)=S(i).yd;%Y������
					cluster=cluster+1;%��ͷ����+1
					 
					 
					distance;
					if (distance>do)%�������ͨ�Ű뾶
						S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); %��������
					end
					if (distance<=do)%����С��ͨ�Ű뾶
						S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Efs*4000*( distance * distance )); %��������
					end
				end 
				
			end
		end 
	end

	 
	 
	 
	STATISTICS(r+1).CLUSTERHEADS=cluster-1;%r�ֺ��ͷ��
	CLUSTERHS(r+1)=cluster-1;%r�ֺ��ͷ��
	 
	 
	for i=1:1:n
		if ( S(i).type=='N' && S(i).E>0 )%ѡ�������ڵ����ش�ͷ
			if(cluster-1>=1)%��ͷ��������2��
				min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );%���ڵ����̾���
				min_dis_cluster=1;%������С�Ĵ�ͷ��
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
				 
				if(min_dis>0)%������ɢ
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
	for i=1:1:n%��ǰ�ڵ���
			if S(i).E > 0%����ڵ�iʣ����������0
				STATISTICS(r+1).ENERGY = STATISTICS(r+1).ENERGY +S(i).E;%r�ֺ�ڵ�ʣ���������Ͻڵ�i��ʣ������
		end
	end
	 
	 
	%[vx,vy]=voronoi(X,Y);
	%plot(X,Y,'r*',vx,vy,'b-');
	%hold on;
	% voronoi(X,Y);
	%axis([0 xm 0 ym]);
	hold off;
end
for i=2:rmax%��ǰ�ڵ���
	mylive(i) = n - STATISTICS(i).DEAD;
	myenergy(i) = STATISTICS(i).ENERGY;%ʣ������
end
mylive(1)=100;
myenergy(1)=S(1).E+(n-1)*Eo;
figure(2);%���ͼ��2
hold on;%��������
plot(mylive,'color','r');%�ú�ɫ������ڵ���
xlabel('������');
ylabel('���ڵ�');
title('���ڵ�ͼ');
figure(3);%���ͼ��3
hold on;%��������
plot(myenergy,'color','r');%�ú�ɫ���ʣ������
xlabel('������');
ylabel('ʣ�������ڵ�');
title('ʣ������ͼ');