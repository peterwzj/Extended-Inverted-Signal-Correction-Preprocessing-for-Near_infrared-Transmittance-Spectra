function [eisc,ai,bi,ci,di,ei] = EISC_2(data,wavelength)
%EISC 扩展反向信号校正算法实现
%输入参数data：原始光谱数据n*p;wavelength:波长：(n,)
%返回值：扩展反向信号校正预处理后的光谱
%syms z w;
[n,p] = size(data);
eisc = zeros(n,p);
m = mean(data,1);%按行塌缩方向求均值
m = m';%1*100 ~ 100*1
w = wavelength;
ai = zeros(1,n);
bi = zeros(1,n);
ci = zeros(1,n);
di = zeros(1,n);
ei = zeros(1,n);
for i = 1:n
    z = data(i,:);
    z2 = z.^2;
    w2 = w'.^2;
    A = [z2;w2;z;w';ones(1,p)];%5*100
    A = A';%100*5
    cc = A\m;%这里使用zuo除运算符，后台自动使用QR分解。
    %原方程是方程数远大于未知数个数（100>5).cc:5*1
    ai(i)=cc(5);
    bi(i)=cc(3);
    ci(i)=cc(1);
    di(i) = cc(4);
    ei(i) = cc(2);%符合文章公式9的排序
    eisc(i,:) = cc(1)*z.^2+cc(2)*w'.^2+cc(3)*z+cc(4)*w'+cc(5);
end