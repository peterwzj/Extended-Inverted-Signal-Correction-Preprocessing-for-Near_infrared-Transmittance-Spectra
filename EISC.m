function eisc = EISC(data,wavelength)
%EISC 扩展反向信号校正算法实现
%输入参数data：原始光谱数据n*p;wavelength:波长：(n,)
%返回值：扩展反向信号校正预处理后的光谱
%syms z w;
[n,p] = size(data);
eisc = zeros(n,p);
m = mean(data,1);%按行塌缩方向求均值
w = wavelength;
for i = 1:n
    z = data(i,:);
    z2 = z.^2;
    w2 = w'.^2;
    A = [z2;w2;z;w';ones(1,p)];
    cc = m/A;
    eisc(i,:) = cc(1)*z.^2+cc(2)*w'.^2+cc(3)*z+cc(4)*w'+cc(5);
end