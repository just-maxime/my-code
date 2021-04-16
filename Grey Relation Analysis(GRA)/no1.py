# -*- coding: utf-8 -*-
"""
Created on Fri May  1 18:42:31 2020

@author: maxim
"""

import pandas as pd
x=pd.read_excel('data.xlsx')
x=x.iloc[:,1:].T

#数据均值化处理
x_mean=x.mean(axis=1)
for i in range(x.index.size):
    x.iloc[i,:] = x.iloc[i,:]/x_mean[i]

#提取参考队列和比较队列
ck=x.iloc[0,:]
cp=x.iloc[1:,:]

# 比较队列与参考队列相减
t=pd.DataFrame()
for j in range(cp.index.size):
    temp=pd.Series(cp.iloc[j,:]-ck)
    t=t.append(temp,ignore_index=True)

#求最大差和最小差
mmax=t.abs().max().max()
mmin=t.abs().min().min()
rho=0.5
#求关联系数
ksi=((mmin+rho*mmax)/(abs(t)+rho*mmax))


#求关联度
r=ksi.sum(axis=1)/ksi.columns.size

#关联度排序
result=r.sort_values(ascending=False)
print(r)