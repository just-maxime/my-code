clc,clear

format long

load('data.txt');%��ԭʼ���ݱ����ڴ��ı��ļ�data.txt��

X=data(:,[1:9]);%XΪ���������Ϊ��������ĸ���

X=X';

Y=data(:,[10:12]);%YΪ���������Ϊ��������ĸ���

Y=Y';           

n=size(X',1);m=size(X,1);s=size(Y,1);

A=[-X' Y'];

b=zeros(n,1);

LB=zeros(m+s,1);
UB=[];

for i=1:n;

  f=[zeros(1,m)  -Y(:,i)'];

  Aeq=[X(:,i)',zeros(1,s)];beq=1;

  w(:,i)=linprog(f,A,b,Aeq,beq,LB,UB);

    E(i,i)=Y(:,i)'*w(m+1:m+s,i);

end

theta=diag(E)';

fprintf('��DEA�����Դ˵�������۽��Ϊ��\n');

disp(theta);

omega=w(1:m,:)

mu=w(m+1:m+s,:)