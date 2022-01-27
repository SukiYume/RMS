import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from sklearn import preprocessing
from sklearn.cluster import DBSCAN

def ArPLS(y, lam=1e4, ratio=0.05, itermax=10):
    N = len(y)
    D = sparse.eye(N, format='csc')
    D = D[1:] - D[:-1]
    D = D[1:] - D[:-1]
    D, w = D.T, np.ones(N)
    for _ in range(itermax):
        W = sparse.diags(w, 0, shape=(N, N))
        Z = W + lam * D.dot(D.T)
        z = spsolve(Z, w * y)
        d = y - z
        dn = d[d < 0]
        m, s = np.mean(dn), np.std(dn)
        wt = 1. / (1 + np.exp(2 * (d - (2 * s - m)) / s))
        if np.linalg.norm(w - wt) / np.linalg.norm(w) < ratio: break
        w = wt
    return z

def get_time_index_new(I, noise_sigma=2, squeeze_frac=3):
    ##### 平滑脉冲轮廓以取最小值
    I_m = ArPLS(I, lam=5, ratio=1, itermax=100)
    ##### 取极小值索引
    extremum_min = (np.diff(np.sign(np.diff(I_m))) > 0).nonzero()[0] + 1
    ##### 计算噪声水平，默认取sigma=2
    max_noise = np.median(I[extremum_min]) + noise_sigma * np.std(I[extremum_min])
    ##### 选择在噪声水平以内的极小值索引，插入最大值索引并排序
    I_index = np.sort([i for i in extremum_min if I[i] < max_noise] + [np.argmax(I_m)])
    ##### 选择最大值索引所在索引
    max_index = np.where(I_index == np.argmax(I_m))[0][0]
    ##### 选择最大值索引左右一个索引作为起终点，并收缩1/5
    if max_index >= 1:
        start, end = I_index[[max_index-1, max_index+1]]
        start, end = start + int((end - start) * 1 / squeeze_frac / 2), end - int((end - start) * 1 / squeeze_frac / 2)
    else:
        start, end = 0, I_index[max_index+1]
    signal_time = np.zeros(len(I))
    signal_time[start: end] = 1
    return signal_time

def get_time_index(I):
    ##### 按时间轴标准化处理数据
    I_scale = preprocessing.MinMaxScaler().fit(I.T).transform(I.T)
    ##### 变换连通半径，对事件通道聚类
    for i in np.arange(0.1, np.linalg.norm(I_scale, axis=1).max()*2, 0.01):
        signal_time = DBSCAN(min_samples=2, eps=i).fit(np.sum(I_scale, axis=1).reshape(-1, 1)).labels_
        if len(set(signal_time)) > 1 and 0<sum(signal_time!=0)<len(signal_time)/3 and sum((np.diff(np.where(signal_time!=0))>2)[0])<1:
            break
    ##### 如果不够长，用新方法取值
    if len(set(signal_time)) <= 2 or len(np.where(signal_time!=0)[0]) < 5:
        signal_time = np.zeros(I.shape[1])
        signal_time[np.argmax(np.sum(I, axis=0))-2: np.argmax(np.sum(I, axis=0))+2] = 1
#        signal_time = get_time_index_new(np.mean(I, axis=0), noise_sigma=2, squeeze_frac=1.5)
    return signal_time

def get_freq_index(I):
    ##### 按频率轴标准化处理数据
    I_scale = preprocessing.MinMaxScaler().fit(I).transform(I)
    ##### 变换连通半径，对频率通道聚类
    for i in np.arange(0.1, np.linalg.norm(I_scale, axis=1).max(), 0.1):
        signal_chan = DBSCAN(min_samples=2, eps=i).fit(I_scale).labels_
        if len(set(signal_chan)) > 1 and 0 < sum(signal_chan!=0) < len(signal_chan) / 2:
            break
    ##### 改变聚类符号，统一之后判断语句
    if sum(signal_chan!=0) < len(signal_chan) / 10:
        signal_chan = -1 * (signal_chan + 1)
        signal_chan[np.sum(I, axis=1)<=10] = 0
    return signal_chan

def get_pulse(I, Q, U, V, freq):
    signal_time = get_time_index(I)
    signal_chan = get_freq_index(I)
    I_clean = I[signal_chan!=0][:, signal_time!=0]
    Q_clean = Q[signal_chan!=0][:, signal_time!=0]
    U_clean = U[signal_chan!=0][:, signal_time!=0]
    V_clean = V[signal_chan!=0][:, signal_time!=0]
    freq_clean = freq[signal_chan!=0]
    off_profile = np.sum(I[:, signal_time==0], axis=0)
    rms = np.random.choice(off_profile, replace=False, size=int(len(off_profile)*0.8)).std()
    snr = np.sum(I_clean) / np.sqrt(I_clean.shape[1]) / rms
    return I_clean, Q_clean, U_clean, V_clean, freq_clean, snr