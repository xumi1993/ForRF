# -*- coding: utf-8 -*-
import numpy as np
import obspy
from obspy.io.sac import SACTrace
from obspy.signal.util import next_pow_2
from math import pi
from scipy.fftpack import fft, ifft, ifftshift
from scipy.signal import fftconvolve, correlate
from scipy.linalg import solve_toeplitz
import matplotlib.pyplot as plt
import time


def gaussFilter(dt, nft, f0):
    df = 1.0 / (nft * dt)
    nft21 = 0.5 * nft + 1
    f = df * np.arange(0, nft21)
    w = 2 * pi * f

    gauss = np.zeros([nft, 1])
    gauss1 = np.exp(-0.25 * (w / f0) ** 2) / dt
    gauss1.shape = (len(gauss1), 1)
    gauss[0:int(nft21)] = gauss1
    gauss[int(nft21):] = np.flipud(gauss[1:int(nft21) - 1])
    gauss = gauss[:, 0]

    return gauss


def gfilter(x, nfft, gauss, dt):
    Xf = fft(x, nfft)
    Xf = Xf * gauss * dt
    xnew = ifft(Xf, nfft).real
    return xnew


def correl(R, W, nfft):
    x = ifft(fft(R, nfft) * np.conj(fft(W, nfft)), nfft)
    x = x.real
    return x


def phaseshift(x, nfft, dt, tshift):
    Xf = fft(x, nfft)

    shift_i = int(tshift / dt)
    p = 2 * pi * np.arange(1, nfft + 1) * shift_i / nfft
    Xf = Xf * np.vectorize(complex)(np.cos(p), -np.sin(p))
    np.savetxt('t2.txt', Xf)
    x = ifft(Xf, nfft) / np.cos(2 * pi * shift_i / nfft)
    x = x.real
    return x


def deconit(uin, win, dt, nt=None, tshift=10, f0=2.0, itmax=400, minderr=0.001, phase='P'):
    """
    Created on Wed Sep 10 14:21:38 2014

    In:
    uin = numerator (radial for PdS)
    win = denominator (vertical component for PdS)
    dt = sample interval (s)
    nt = number of samples
    tshift = Time until beginning of receiver function (s)
    f0 = width of gaussian filter
    itmax = max # iterations
    minderr = Min change in error required for stopping iterations

    Out:
    RFI = receiver function
    rms = Root mean square error for predicting numerator after each iteration

    @author: Mijian Xu @ NJU
    """
    # print('Iterative Decon (Ligorria & Ammon):\n')
    if len(uin) != len(win):
        raise ValueError('The two input trace must be in same length')
    elif nt is None:
        nt = len(uin)
    else:
        pass

    rms = np.zeros(itmax)
    nfft = next_pow_2(nt)
    p0 = np.zeros(nfft)

    u0 = np.zeros(nfft)
    w0 = np.zeros(nfft)

    u0[0:nt] = uin
    w0[0:nt] = win

    gaussF = gaussFilter(dt, nfft, f0)
    # gaussF = _gauss_filter(dt, nfft, f0)

    u_flt = gfilter(u0, nfft, gaussF, dt)
    w_flt = gfilter(w0, nfft, gaussF, dt)


    wf = fft(w0, nfft)
    r_flt = u_flt

    powerU = np.sum(u_flt ** 2)

    it = 0
    sumsq_i = 1
    d_error = 100 * powerU + minderr
    maxlag = 0.5 * nfft
    # print('\tMax Spike Display is ' + str((maxlag) * dt))
    while np.abs(d_error) > minderr and it < itmax:
        rw = correl(r_flt, w_flt, nfft)
        rw = rw / np.sum(w_flt ** 2)

        if phase == 'P':
            i1 = np.argmax(np.abs(rw[0:int(maxlag) - 1]))
        else:
            i1 = np.argmax(np.abs(rw))
        amp = rw[i1] / dt

        # print(i1)

        p0[i1] = p0[i1] + amp
        p_flt = gfilter(p0, nfft, gaussF, dt)
        p_flt = gfilter(p_flt, nfft, wf, dt)

        r_flt = u_flt - p_flt
        sumsq = np.sum(r_flt ** 2) / powerU
        rms[it] = sumsq
        d_error = 100 * (sumsq_i - sumsq)

        sumsq_i = sumsq

        it = it + 1
    p_flt = gfilter(p0, nfft, gaussF, dt)
    p_flt = phaseshift(p_flt, nfft, dt, tshift)

    RFI = p_flt[0:nt]
    rms = rms[0:it - 1]

    return RFI, rms, it 


def dorf(r, z):
    # nft = next_pow_2(r.npts)
    # ga = gaussFilter(r.delta, nft, 2.0)
    # l = np.flip(ldata.data, axis=0)
    # q = np.flip(qdata.data, axis=0)
    time_axis = np.arange(r.npts)*r.delta+r.b
    time_start=time.time()
    rf, rms, it = deconit(r.data, z.data, r.delta, tshift=-r.b, f0=2, itmax=200)
    time_end=time.time()
    print('time cost',time_end-time_start,'s')
    plt.plot(time_axis, rf, lw=3, label='seispy')
    rr = SACTrace.read('rf.sac')
    plt.plot(time_axis, rr.data, '--', lw=1.5, label='ForRF')
    plt.legend()
    plt.xlim(-2,30)
    plt.xlabel('Time (s)')
    plt.savefig('cmp.png', dpi=300)


def dofilter(r, f0=2):
    nft = next_pow_2(r.npts)
    ga = gaussFilter(r.delta, nft, f0)
    data = gfilter(r.data, nft, ga, r.delta)
    np.savetxt('pyflt.txt',data)


def dofft(r):
    nft = next_pow_2(r.npts)
    z = fft(r.data, nft)
    z = ifft(z, nft)
    np.savetxt('pyfft.txt', z)

if __name__ == '__main__':
    r = SACTrace.read('test_R.sac')
    z = SACTrace.read('test_Z.sac')

    # dofilter(r)
    dorf(r,z)
    # dofft(r)
