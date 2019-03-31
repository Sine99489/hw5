第五次作业Matlab代码

~~~matlab
function g=im_dftfilt(f,H)
%图像频域滤波函数
f=im2double(f);
[m,n]=size(f);
F=fft2(f,2*m,2*n);
g=ifft2(H.*F);
g=g(1:m,1:n);
g=im2uint8(g);
~~~

~~~matlab
function D=getD(m,n)
%求距离矩阵
is1=0:m-1;
is2=0:n-1;
i1=find(is1>m/2);
i2=find(is2>n/2);
is1(i1)=is1(i1)-m;
is2(i2)=is2(i2)-n;
[y,x]=meshgrid(is2,is1);
D=hypot(x,y);
~~~

~~~matlab
function H=getfilter(type,s1,s2,D0,n)
%获取滤波器函数
D=getD(s1,s2);
switch type
    case 'Gaussian'
        H=exp(-(D.^2)./(2*(D0^2)));
    case 'Butterworth'
        H=1./(1+(D/D0).^(2*n));
    case 'Laplace'
        H=-4*(pi^2)*(D.^2);
       % H=H/H(floor(s1/2),floor(s2/2));
    otherwise
        error('Unknow type!')
end
~~~

~~~matlab
%p5_1.m
%频域低通滤波器
%设计低通滤波器包括 butterworth and Gaussian (选择合适的半径，计算功率谱比),平滑测试图像test1和2
test1=imread('test1.pgm');
test2=imread('test2.tif');
figure(1)
subplot(2,4,5)
imshow(test1);
title('test1原图像');
figure(2)
subplot(2,4,5)
imshow(test2);
title('test2原图像');
figure(3)
subplot(2,4,5)
imshow(test1);
title('test1原图像');
figure(4)
subplot(2,4,5)
imshow(test2);
title('test2原图像');

s1=size(test1);
s2=size(test2);
D0_1=round(s1(1)*[5,10,20]/100);
D0_2=round(s2(1)*[5,10,20]/100);

%图像离散傅里叶变换
test1=im2double(test1);
test2=im2double(test2);
T1=fft2(test1,2*s1(1),2*s1(2));
T2=fft2(test2,2*s2(1),2*s2(2));

%计算总功率
P1=sum(abs(T1(:)).^2);
P2=sum(abs(T2(:)).^2);

%滤波处理
for i=1:3
    %butterworth滤波处理
    Hb1=getfilter('Butterworth',2*s1(1),2*s1(2),D0_1(i),2);
    Hb2=getfilter('Butterworth',2*s2(1),2*s2(2),D0_2(i),2);
    Gb1=Hb1.*T1;
    gb1=ifft2(Gb1);
    gb1=gb1(1:s1(1),1:s1(2));
    Gb2=Hb2.*T2;
    gb2=ifft2(Gb2);
    gb2=gb2(1:s2(1),1:s2(2));
    
    %高斯滤波处理
    Hg1=getfilter('Gaussian',2*s1(1),2*s1(2),D0_1(i));
    Hg2=getfilter('Gaussian',2*s2(1),2*s2(2),D0_2(i));
    Gg1=Hg1.*T1;
    gg1=ifft2(Gg1);
    gg1=gg1(1:s1(1),1:s1(2));
    Gg2=Hg2.*T2;
    gg2=ifft2(Gg2);
    gg2=gg2(1:s2(1),1:s2(2));
    
    %计算功率谱比
    rb1=roundn(sum(abs(Gb1(:)).^2)/P1*100,-2);
    rb2=roundn(sum(abs(Gb2(:)).^2)/P2*100,-2);
    rg1=roundn(sum(abs(Gg1(:)).^2)/P1*100,-2);
    rg2=roundn(sum(abs(Gg2(:)).^2)/P2*100,-2);
    
    %显示图像
    figure(1);
    subplot(2,4,i+1)
    imshow(fftshift(Hb1),[]);
    title(['Butterworth,D0=',num2str(D0_1(i))]);
    subplot(2,4,i+5)
    imshow(gb1);
    title(['功率谱比为',num2str(rb1),'%'])
    
    figure(2);
    subplot(2,4,i+1)
    imshow(fftshift(Hb2),[]);
    title(['Butterworth,D0=',num2str(D0_2(i))]);
    subplot(2,4,i+5)
    imshow(gb2);
    title(['功率谱比为',num2str(rb2),'%']);
    
    figure(3);
    subplot(2,4,i+1)
    imshow(fftshift(Hg1),[]);
    title(['Gaussian,D0=',num2str(D0_1(i))]);
    subplot(2,4,i+5)
    imshow(gg1);
    title(['功率谱比为',num2str(rg1),'%']);
    
    figure(4);
    subplot(2,4,i+1)
    imshow(fftshift(Hg2),[]);
    title(['Gaussian,D0=',num2str(D0_2(i))]);
    subplot(2,4,i+5)
    imshow(gg2);
    title(['功率谱比为',num2str(rg2),'%']);
end
~~~

~~~matlab
%p5_2.m
%频域高通滤波器
%设计高通滤波器包括butterworth and Gaussian，在频域增强边缘。选择半径和计算功率谱比，测试图像test3,4
test3=imread('test3_corrupt.pgm');
test4=imread('test4 copy.bmp');
figure(1)
subplot(2,4,5)
imshow(test3);
title('test3原图像');
figure(2)
subplot(2,4,5)
imshow(test4);
title('test4原图像');
figure(3)
subplot(2,4,5)
imshow(test3);
title('test3原图像');
figure(4)
subplot(2,4,5)
imshow(test4);
title('test4原图像');
s1=size(test3);
s2=size(test4);
D0_1=round(s1(1)*[5,10,20]/100);
D0_2=round(s2(1)*[5,10,20]/100);

%图像离散傅里叶变换
test3=im2double(test3);
test4=im2double(test4);
T1=fft2(test3,2*s1(1),2*s1(2));
T2=fft2(test4,2*s2(1),2*s2(2));

%计算总功率
P1=sum(abs(T1(:)).^2);
P2=sum(abs(T2(:)).^2);

%滤波处理
for i=1:3
    %butterworth滤波处理
    Hb1=1-getfilter('Butterworth',2*s1(1),2*s1(2),D0_1(i),2);
    Hb2=1-getfilter('Butterworth',2*s2(1),2*s2(2),D0_2(i),2);
    Gb1=Hb1.*T1;
    gb1=ifft2(Gb1);
    gb1=gb1(1:s1(1),1:s1(2));
    Gb2=Hb2.*T2;
    gb2=ifft2(Gb2);
    gb2=gb2(1:s2(1),1:s2(2));
    
    %高斯滤波处理
    Hg1=1-getfilter('Gaussian',2*s1(1),2*s1(2),D0_1(i));
    Hg2=1-getfilter('Gaussian',2*s2(1),2*s2(2),D0_2(i));
    Gg1=Hg1.*T1;
    gg1=ifft2(Gg1);
    gg1=gg1(1:s1(1),1:s1(2));
    Gg2=Hg2.*T2;
    gg2=ifft2(Gg2);
    gg2=gg2(1:s2(1),1:s2(2));
    
    %计算功率谱比
    rb1=roundn(sum(abs(Gb1(:)).^2)/P1*100,-2);
    rb2=roundn(sum(abs(Gb2(:)).^2)/P2*100,-2);
    rg1=roundn(sum(abs(Gg1(:)).^2)/P1*100,-2);
    rg2=roundn(sum(abs(Gg2(:)).^2)/P2*100,-2);
    
     %显示图像
    figure(1);
    subplot(2,4,i+1)
    imshow(fftshift(Hb1),[]);
    title(['Butterworth,D0=',num2str(D0_1(i))]);
    subplot(2,4,i+5)
    imshow(gb1,[]);
    title(['功率谱比为',num2str(rb1),'%'])
    
    figure(2);
    subplot(2,4,i+1)
    imshow(fftshift(Hb2),[]);
    title(['Butterworth,D0=',num2str(D0_2(i))]);
    subplot(2,4,i+5)
    imshow(gb2,[]);
    title(['功率谱比为',num2str(rb2),'%']);
    
    figure(3);
    subplot(2,4,i+1)
    imshow(fftshift(Hg1),[]);
    title(['Gaussian,D0=',num2str(D0_1(i))]);
    subplot(2,4,i+5)
    imshow(gg1,[]);
    title(['功率谱比为',num2str(rg1),'%']);
    
    figure(4);
    subplot(2,4,i+1)
    imshow(fftshift(Hg2),[]);
    title(['Gaussian,D0=',num2str(D0_2(i))]);
    subplot(2,4,i+5)
    imshow(gg2,[]);
    title(['功率谱比为',num2str(rg2),'%']);
end
~~~

~~~matlab
%p5_3.m
%其他高通滤波器
%拉普拉斯和Unmask，对测试图像test3,4滤波
test3=imread('test3_corrupt.pgm');
test4=imread('test4 copy.bmp');
s1=size(test3);
s2=size(test4);

%图像离散傅里叶变换
test3=im2double(test3);
test4=im2double(test4);
T1=fft2(test3,2*s1(1),2*s1(2));
T2=fft2(test4,2*s2(1),2*s2(2));

%拉普拉斯算子锐化处理
Hl1=getfilter('Laplace',2*s1(1),2*s1(2));
Hl2=getfilter('Laplace',2*s2(1),2*s2(2));
gl1=ifft2(Hl1.*T1);
gl1=gl1(1:s1(1),1:s1(2));
gl1=test3/max(test3(:))-gl1/max(abs(gl1(:)));
gl2=ifft2(Hl2.*T2);
gl2=gl2(1:s2(1),1:s2(2));
gl2=test4/max(test4(:))-gl2/max(abs(gl2(:)));

%非锐化掩蔽处理
Hg1=getfilter('Gaussian',2*s1(1),2*s1(2),round(0.05*s1(1)*2));
Hg2=getfilter('Gaussian',2*s2(1),2*s2(2),round(0.05*s2(1)*2));
gg1=ifft2(Hg1.*T1);
gu1=2*test3-gg1(1:s1(1),1:s1(2));
gg2=ifft2(Hg2.*T2);
gu2=2*test4-gg2(1:s2(1),1:s2(2));

%显示图像
figure(1);
subplot(1,3,1)
imshow(test3);
title(['test3原图像']);
subplot(1,3,2)
imshow(gl1);
title(['用拉普拉斯算子锐化test3图像']);
subplot(1,3,3)
imshow(gu1);
title(['非锐化掩蔽处理test3图像']);

figure(2);
subplot(1,3,1)
imshow(test4);
title(['test4原图像']);
subplot(1,3,2)
imshow(gl2);
title(['用拉普拉斯算子锐化test4图像']);
subplot(1,3,3)
imshow(gu2);
title(['非锐化掩蔽处理test4图像']);

~~~

