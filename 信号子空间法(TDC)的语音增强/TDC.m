%基于Ephraim和Van Trees提出的信号子空间法(TDC)的语音增强程序，适用于噪声为白噪声的情况
clear all;
%--------------------------------参数定义----------------------------------
frame_len=80;              %帧长
step_len=0.5*frame_len;    %分帧时的步长，相当于重叠50%
N=2;                       %计算Toeplitz协方差矩阵时用到的前后相邻的帧数,N为偶数
v=0.05;                    %噪声抑制系数，推荐v=2或3;
u=0.05;
fs=44100;
%-------------------------------读入带噪语音文件----------------------------
[tidy,nbits]=audioread('before.wav');%纯净语音文件
[wavin,nbits]=audioread('add_noise.wav');%加噪语音文件
wav_length=length(wavin);
%--------------------------------分帧--------------------------------------
R = step_len;
L = frame_len; 
f = (wav_length-mod(wav_length,frame_len))/frame_len;
k = 2*f-1;                          % 帧数        
frame_num=k;
for r = 1:k 
    y = wavin(1+(r-1)*R:L+(r-1)*R); % 对带噪语音帧间重叠一半取值；
    out(1+(r-1)*L:r*L) = y(1:L);    % 得到一个新帧数的序列
end
inframe=reshape(out,frame_len,k);   % 改变序列的形状
%-----------------------------子空间法语音增强------------------------------
%求噪声的Toeplitz协方差矩阵，Rn=n_var*I;
n_var=var(wavin(1:2200),1);  %取前3000采样点用于估计噪声方差
xv=n_var;                    %定义Rx的特征值判别阈值
Rn=n_var*eye(frame_len);
Ry=zeros(frame_len,frame_len);          %定义帧信号的Toeplitz协方差矩阵
seq=zeros(1,frame_len);                 %定义相邻6帧的自相关序列，长度等于帧长
L=N*frame_len;                          %定义相邻N帧的长度
outframe=zeros(frame_len,frame_num);    %定义增强后的矩阵
%对每一帧进行处理，得到增强后的信号,因为计算Toeplitz协方差矩阵要用到N帧，其中当前
%帧处在第(N/2+1)帧（注意：计算Toeplitz协方差矩阵的帧是不重叠的，而inframe相邻帧
%之间重叠50%），因此，inframe的前N帧和最后(N/2-1)*2帧无法估计Toeplitz协方差矩阵，
%且这(N-1)*2帧一般是噪声信号，故不必对其处理，直接将out_frame的这些帧置0
%------------------------估计噪声的特征值----------------------------------
%  Ann=noise_evaluate_min(wavin,frame_len,frame_num,N,step_len,xv);
for i=(N+1):(frame_num-N-2)
%1.构造当前帧的Toeplitz协方差矩阵
    for j=0:(frame_len-1)
        bgn_point=(i-N-1)*step_len+1;       %相邻6帧的起始点
        end_point=(i+N-1)*step_len;         %相邻6帧的终点
        seq(j+1)=wavin(bgn_point:(end_point-j))'*...
                 wavin((bgn_point+j):end_point)/L;
    end;
    Ry=toeplitz(seq);
%2.对Ry进行特征分解，求它的特征向量矩阵Uy和特征值矩阵Ay
    [Uy,Ay]=eig(Ry);% [Uy,Ay]=eig(Ry,'nobalance');
%3.求出Ay中大于阈值xv的个数M，即为信号子空间维数，并由此求出语音子空间的特征值矩
%阵Ax=Ay-Rn(维数:M*M),以及对应的特征向量矩阵Ux(维数：frame_len*M)
    [I,J]=find(Ay>xv);
    M=length(I);
    Ax=zeros(M,M);
    Ux=zeros(frame_len,M);
    Ay_seq=zeros(1,frame_len);
    A=sort(Ay);        
    seq1=A(frame_len,:);
    [Ay_seq,IX]=sort(seq1);
%------------------------根据子空间能量公式计算winner比值G_k-----------------   
%       e1(i)=sum(Ay_seq(frame_len-M+1:frame_len));        % 信号能量
%       e2(i)=frame_len*n_var;                             % 噪声能量
%       G_k(i)=e1(i)/(e1(i)+e2(i));   
%----------------------winner----------------------------------------------
     for k=1:M
        num=frame_len-k+1;
%             S1=Ay_seq(num)-Ann(num);
%             S=max(S1,0);
%             Ax=G_k(i).*S;
%             Ax=Ay_seq(num)-Ann(num);
%             Ax(k,k)=G_k(i).*S(num);
              Ax(k,k)=Ay_seq(num)-xv;
              Ux(:,k)=Uy(:,IX(num));
    end;
    Gu=Ax./(Ax+u*n_var);
    outframe(:,i)=Ux*Gu*(conj(Ux)')*inframe(:,i);
end;
%----------------------将增强后的帧信号连接成语音---------------------------
wavout=zeros(1,(frame_num-1)*step_len+frame_len);
for t=1:frame_num
    num1=(t-1)*step_len+1;
    num2=(t+1)*step_len;
    wavout(num1:num2)=wavout(num1:num2)+(hamming(frame_len).*outframe(:,t))';
end;
wavout=wavout';
audiowrite('after.wav',wavout,44100);
%-----------------------信噪比--------------------------------------------
x100=audioread('before.wav');
y100=audioread('add_noise.wav');
p100=norm(x100).^2;
p200=norm(x100-y100).^2;
SNR_before=10*log(p100/p200);

x100=audioread('before.wav');
y200=audioread('after.wav');
p20=norm(x100-y200).^2;
SNR_after=10*log(p100/p20);


%audiowrite(wavout,fs,nbits,['EVT_' num2str(frame_len) '_' num2str(N) '_51u_' filename]);
%-----------------------将处理前后的结果进行作图比较-------------------------
figure(1);
subplot(3,1,1);plot(tidy);xlabel('(a)原始语音（采样点数）');ylabel('幅度'); axis([2.5*10^4 2.8*10^4 -0.03 0.03]);
subplot(3,1,2);plot(wavin);xlabel('(b)带噪语音（采样点数）');ylabel('幅度');axis([2.5*10^4 2.8*10^4 -0.03 0.03]);
subplot(3,1,3);plot(wavout);xlabel('(c)子空间法增强语音-TDC（采样点数）');ylabel('幅度');axis([2.5*10^4 2.8*10^4 -0.03 0.03]);
