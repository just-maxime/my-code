%����Ephraim��Van Trees������ź��ӿռ䷨(TDC)��������ǿ��������������Ϊ�����������
clear all;
%--------------------------------��������----------------------------------
frame_len=80;              %֡��
step_len=0.5*frame_len;    %��֡ʱ�Ĳ������൱���ص�50%
N=2;                       %����ToeplitzЭ�������ʱ�õ���ǰ�����ڵ�֡��,NΪż��
v=0.05;                    %��������ϵ�����Ƽ�v=2��3;
u=0.05;
fs=44100;
%-------------------------------������������ļ�----------------------------
[tidy,nbits]=audioread('before.wav');%���������ļ�
[wavin,nbits]=audioread('add_noise.wav');%���������ļ�
wav_length=length(wavin);
%--------------------------------��֡--------------------------------------
R = step_len;
L = frame_len; 
f = (wav_length-mod(wav_length,frame_len))/frame_len;
k = 2*f-1;                          % ֡��        
frame_num=k;
for r = 1:k 
    y = wavin(1+(r-1)*R:L+(r-1)*R); % �Դ�������֡���ص�һ��ȡֵ��
    out(1+(r-1)*L:r*L) = y(1:L);    % �õ�һ����֡��������
end
inframe=reshape(out,frame_len,k);   % �ı����е���״
%-----------------------------�ӿռ䷨������ǿ------------------------------
%��������ToeplitzЭ�������Rn=n_var*I;
n_var=var(wavin(1:2200),1);  %ȡǰ3000���������ڹ�����������
xv=n_var;                    %����Rx������ֵ�б���ֵ
Rn=n_var*eye(frame_len);
Ry=zeros(frame_len,frame_len);          %����֡�źŵ�ToeplitzЭ�������
seq=zeros(1,frame_len);                 %��������6֡����������У����ȵ���֡��
L=N*frame_len;                          %��������N֡�ĳ���
outframe=zeros(frame_len,frame_num);    %������ǿ��ľ���
%��ÿһ֡���д����õ���ǿ����ź�,��Ϊ����ToeplitzЭ�������Ҫ�õ�N֡�����е�ǰ
%֡���ڵ�(N/2+1)֡��ע�⣺����ToeplitzЭ��������֡�ǲ��ص��ģ���inframe����֡
%֮���ص�50%������ˣ�inframe��ǰN֡�����(N/2-1)*2֡�޷�����ToeplitzЭ�������
%����(N-1)*2֡һ���������źţ��ʲ��ض��䴦��ֱ�ӽ�out_frame����Щ֡��0
%------------------------��������������ֵ----------------------------------
%  Ann=noise_evaluate_min(wavin,frame_len,frame_num,N,step_len,xv);
for i=(N+1):(frame_num-N-2)
%1.���쵱ǰ֡��ToeplitzЭ�������
    for j=0:(frame_len-1)
        bgn_point=(i-N-1)*step_len+1;       %����6֡����ʼ��
        end_point=(i+N-1)*step_len;         %����6֡���յ�
        seq(j+1)=wavin(bgn_point:(end_point-j))'*...
                 wavin((bgn_point+j):end_point)/L;
    end;
    Ry=toeplitz(seq);
%2.��Ry���������ֽ⣬������������������Uy������ֵ����Ay
    [Uy,Ay]=eig(Ry);% [Uy,Ay]=eig(Ry,'nobalance');
%3.���Ay�д�����ֵxv�ĸ���M����Ϊ�ź��ӿռ�ά�������ɴ���������ӿռ������ֵ��
%��Ax=Ay-Rn(ά��:M*M),�Լ���Ӧ��������������Ux(ά����frame_len*M)
    [I,J]=find(Ay>xv);
    M=length(I);
    Ax=zeros(M,M);
    Ux=zeros(frame_len,M);
    Ay_seq=zeros(1,frame_len);
    A=sort(Ay);        
    seq1=A(frame_len,:);
    [Ay_seq,IX]=sort(seq1);
%------------------------�����ӿռ�������ʽ����winner��ֵG_k-----------------   
%       e1(i)=sum(Ay_seq(frame_len-M+1:frame_len));        % �ź�����
%       e2(i)=frame_len*n_var;                             % ��������
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
%----------------------����ǿ���֡�ź����ӳ�����---------------------------
wavout=zeros(1,(frame_num-1)*step_len+frame_len);
for t=1:frame_num
    num1=(t-1)*step_len+1;
    num2=(t+1)*step_len;
    wavout(num1:num2)=wavout(num1:num2)+(hamming(frame_len).*outframe(:,t))';
end;
wavout=wavout';
audiowrite('after.wav',wavout,44100);
%-----------------------�����--------------------------------------------
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
%-----------------------������ǰ��Ľ��������ͼ�Ƚ�-------------------------
figure(1);
subplot(3,1,1);plot(tidy);xlabel('(a)ԭʼ����������������');ylabel('����'); axis([2.5*10^4 2.8*10^4 -0.03 0.03]);
subplot(3,1,2);plot(wavin);xlabel('(b)��������������������');ylabel('����');axis([2.5*10^4 2.8*10^4 -0.03 0.03]);
subplot(3,1,3);plot(wavout);xlabel('(c)�ӿռ䷨��ǿ����-TDC������������');ylabel('����');axis([2.5*10^4 2.8*10^4 -0.03 0.03]);
