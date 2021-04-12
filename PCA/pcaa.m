clc,clear
 x=xlsread('pca.xlsx','A2:L20');
%x=x';
stdr = std(x);%求各变量标准差
[n,m]=size(x);
sddata = x./stdr(ones(n,1),:);%标准化

[p,princ,egenvalue] = pca(sddata);%调用主成分分析程序
 p3=p(:,1:3)%输出前3个主成分系数，修改，matlab从1开始计数
 sc=princ(:,1:3)%输出3个主成分得分，修改，matlab从1开始计数
egenvalue %输出特征根
per = 100* egenvalue/sum(egenvalue)%输出各个主成分贡献率
