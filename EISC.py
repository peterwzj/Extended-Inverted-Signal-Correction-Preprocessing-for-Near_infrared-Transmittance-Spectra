import numpy as np
#反向信号校正：
def ISC(data):#data为原始光谱
    n,p = data.shape
    isc = np.zeros((n,p))
    mean = np.mean(data,axis=0)
    for i in range(n):
        z = data[i,:]
        coefficients = np.polyfit(z,mean,2)
        isc[i,:] = coefficients[0]*z**2+coefficients[1]*z+coefficients[2]
    return isc

#重写EISC算法函数：使用np.linalg.lstsq来计算超定方程组的解（未知数个数远少于方程个数）
def EISC(data,wavelength):
    n,p = data.shape # data:为原始光谱
    eisc = np.zeros((n,p))
    m = np.mean(data,axis=0)#列方向塌缩求光谱均值,1*100
    m = m.T;#100*1
    w = wavelength.T;#100*1~1*100
    ai = np.zeros((n,))
    bi = np.zeros((n,))
    ci = np.zeros((n,))
    di = np.zeros((n,))
    ei = np.zeros((n,))#文章公式9的系数
    for i in range(n):
        z = data[i,:]
        z2 = z**2
        w2 = w**2
        A = np.vstack([z2,w2,z,w,np.ones((1,p))]) #5*100
        A = A.T #100*5
        cc,residuals,rank,s = np.linalg.lstsq(A,m,rcond=None)
        ai[i] = cc[4];
        bi[i] = cc[2];
        ci[i] = cc[0];
        di[i] = cc[3];
        ei[i] = cc[1];
        eisc[i,:] = cc[0]*z2+cc[1]*w2+cc[2]*z+cc[3]*w+cc[4]
    return eisc,ai,bi,ci,di,ei
