clc,clear
 x=xlsread('pca.xlsx','A2:L20');
%x=x';
stdr = std(x);%���������׼��
[n,m]=size(x);
sddata = x./stdr(ones(n,1),:);%��׼��

[p,princ,egenvalue] = pca(sddata);%�������ɷַ�������
 p3=p(:,1:3)%���ǰ3�����ɷ�ϵ�����޸ģ�matlab��1��ʼ����
 sc=princ(:,1:3)%���3�����ɷֵ÷֣��޸ģ�matlab��1��ʼ����
egenvalue %���������
per = 100* egenvalue/sum(egenvalue)%����������ɷֹ�����
