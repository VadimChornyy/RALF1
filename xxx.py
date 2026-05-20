#https://www.academia.edu/124716825/SOFTWARE_CODE_AND_MATHEMATICAL_ALGORITHM_TO_PREDICT_FUTURE_ATTACKS_VADIM_CHERNOV_KYIV_UKRAINE
import pylab as plt #CIA Submission Reference ID EJJ3C1V1,XULK9Y6R, 1SNGWLU1
import numpy as np #Patent 2011612714
import urllib.request, json 
import pandas as pd
from PIL import Image  
import cv2 as cv
import time as tm 
import hickle as hkl
import os
import multiprocessing as mp
import warnings
from scipy import stats as scp
import dateutil.parser
from operator import itemgetter
import dill 
import csv
from random import sample 
from scipy.signal import savgol_filter
from random import Random
import time
import struct
import scipy.interpolate as interp
import win32api,win32process,win32con         
priorityclasses = [win32process.IDLE_PRIORITY_CLASS,
                win32process.BELOW_NORMAL_PRIORITY_CLASS,
                win32process.NORMAL_PRIORITY_CLASS,
                win32process.ABOVE_NORMAL_PRIORITY_CLASS,
                win32process.HIGH_PRIORITY_CLASS,
                win32process.REALTIME_PRIORITY_CLASS]  

wrkdir = r""
GetSCV=0
Numproc=mp.cpu_count()-1
MxTime=0.66*60*60 # 1/2 haurs

class BaseRandom( Random ):

    def __init__(self, _seed=None):

        super().__init__( _seed )  ## this call creates attribute self._value and sets it

    def seed(self, _seed=None):

        try:
            self.setstate( _seed )
        except:
            super().seed( _seed )

    def __call__(self, _max=1.0):

        return self.uniform( 0.0, _max )
    
class BaseLCG( BaseRandom ):

    def __init__(self, _seedState=None):

        super().__init__( _seedState ) # this call creates attribute self._value and sets it

    def random(self):

        raise NotImplementedError

    def getstate(self):

        return self._value

    def setstate(self, _state):

        raise NotImplementedError

class FastRand63( BaseLCG ):
#Copyright (c) 2016-2019 Philippe Schmouker, schmouk (at) typee.ovh
    def random(self):

        self._value = (9219741426499971445 * self._value + 1) & 0x7fffffffffffffff
        return self._value / 9223372036854775808.0

    def setstate(self, _state):

        if isinstance( _state, int ):
            self._value = _state & 0x7fffffffffffffff
            
        elif isinstance( _state, float ):
            # transforms passed initial seed from float to integer
            if _state < 0.0 :
                _state = -_state
            
            if _state >= 1.0:
                self._value = int(_state+0.5) & 0x7fffffffffffffff
            else:
                self._value = int(_state*0x8000000000000000) & 0x7fffffffffffffff
        
        else:
            t = int(time.time() * 1000.0)
            self._value = t & 0xffffffff
            self._value += (t & 0xff000000) <<  8
            self._value += (t & 0x00ff0000) << 24
            self._value += (t & 0x0000ff00) << 40
            self._value += (t & 0x000000fe) << 63

class BaseLFib64( BaseRandom ):

    def __init__(self, _seedState=None):

        super().__init__( _seedState )
        
    def random(self):

        raise NotImplementedError
            
    def getstate(self):
 
        return (self._list[:], self._index)
            
    def setstate(self, _seedState):

        if not isinstance( _seedState, tuple ):
            self._initIndex( 0 )
            self._initList( _seedState )
            
        elif len( _seedState ) < 2:
            self._initIndex( 0 )
            self._initList( _seedState[0] )
            
        else:
            self._initIndex( _seedState[1] )
            self._list = _seedState[0][:]
                       
    def _initIndex(self, _index):

        try:
            self._index = int( _index ) % self._LIST_SIZE
        except:
            self._index = 0
                       
    def _initList(self, _initialSeed=None):

        myRand = FastRand63( _initialSeed )

        def _getValue( _dummy ):
            myRand()
            v = myRand._value << 1
            return v + (1 if myRand() >= 0.5 else 0)

        self._list = list( map( _getValue, range(self._LIST_SIZE) ) )

class LFib1340( BaseLFib64 ):
    #Copyright (c) 2016-2019 Philippe Schmouker, schmouk (at) typee.ovh

    _LIST_SIZE = 1279 # this 'LFib(2^64, 1279, 861, +)' generator is based on a suite containing 1279 integers

    def random(self):
        """
        This is the core of the pseudo-random generator.
        Returned values are within [0.0, 1.0).
        """
        # evaluates indexes in suite for the i-861 and i-1279 -th values
        k861 = self._index-861
        if k861 < 0:
            k861 += LFib1340._LIST_SIZE
        
        # then evaluates current value
        myValue = (self._list[k861] + self._list[self._index]) & 0xffffffffffffffff
        self._list[self._index] = myValue
        
        # next index
        self._index = (self._index+1) % LFib1340._LIST_SIZE
        
        # then returns float value within [0.0, 1.0)
        return  myValue / 18446744073709551616.0

def XFilterB(dQ2,key=1,key2=1,key3=1):
    #Copyright (c) 2011 Vadim Chernov, ptnt. 2011612714RF
    if key>0:     
        Np=len(dQ2)
        Nf=len(dQ2[0])
        SdQj=np.ones(Np,float)   
        znakSdQj=np.ones(Np,float) 
        
        if key2>0:
            srdQ2=np.mean(dQ2,axis=0)
            dQ2=dQ2-srdQ2
            for i in range(Np):
                znakSdQj[i]=np.mean(dQ2[i])        
            znakSdQj=((znakSdQj<0)-0.5)*2*(-1)
            
        SSdQ=np.std(dQ2)
        
        for i in range(Np):
            if key3>0:
                SdQj[i]=np.std(dQ2[i])
            dQ2[i] = dQ2[i] * znakSdQj[i]            
            if SdQj[i]==0:
                dQ2[i]=np.zeros(Nf,float)
            else:
                dQ2[i]=dQ2[i] / SdQj[i]
   
        SdQ=np.zeros(Nf,float)         
        for j in range(Nf):
            SdQ[j]=np.mean(dQ2[:,j])  
        sSdQ=np.std(SdQ)
        for i in range(Np):
            SdQj_ = np.std(dQ2[i] - SdQ)
            SdQj__ = np.std(dQ2[i])            
            if SdQj__ >0 and sSdQ>0:
                dQ2[i] = dQ2[i] +SdQ * (SdQj_ / sSdQ - 1)
                # dQ2[i] = dQ2[i] * (1+sSdQ / SdQj__) / 2
                dQ2[i] = dQ2[i] *znakSdQj[i] #* SdQj[i]
            else:
                dQ2[i]=np.zeros(Nf,float)    
            
        for i in range(Np):
            dQ2[i] = dQ2[i]* SdQj[i]  
             
        SSdQ_=np.std(dQ2)
        dQ2=dQ2*SSdQ/SSdQ_ 
        
        if key2>0:
            for i in range(Np):
                dQ2[i]=dQ2[i]+srdQ2    
                                
    return dQ2

MaxTemp=90

t_cache = {}  # Global cache for T-Matrix and Psi

import numpy as np
import logging
from scipy.signal import welch
from scipy.stats import wasserstein_distance
from sklearn.covariance import LedoitWolf
import pywt

# === CONFIGURATION ===
USE_KOLMOGOROV_ENHANCEMENTS = True
FAST_MODE = True
KOLMOGOROV_WEIGHTS = {
    'kc': 0.15, 'spectral': 0.05, 'wavelet': 0.10, 'mahal': 0.15,
    'kan': 0.15, 'graph': 0.10, 'wass': 0.05, 'ae': 0.10, 'hurst': 0.05
}

# === VECTORIZED & FAST VERSIONS ===
def fast_lz_complexity(dd, bins=5):
    n, m = dd.shape
    if m < 4: return np.zeros(n)
    mins = dd.min(axis=1, keepdims=True)
    maxs = dd.max(axis=1, keepdims=True) + 1e-12
    binned = ((dd - mins) / (maxs - mins) * (bins-1)).astype(int)
    complexity = np.zeros(n)
    for i in range(n):
        row = binned[i]
        if m < 3:
            complexity[i] = 1.0
            continue
        grams = np.lib.stride_tricks.sliding_window_view(row, 3)
        complexity[i] = len(np.unique(grams, axis=0))
    return np.clip(complexity / (m//3 + 1), 0, 1)

def fast_spectral_entropy(dd, fs=1.0):
    n, m = dd.shape
    if m < 8: return np.zeros(n)
    entropy = np.zeros(n)
    segment = min(64, m//4)
    for i in range(n):
        f, Pxx = welch(dd[i], fs=fs, nperseg=segment, noverlap=0)
        Pxx = Pxx[Pxx > 1e-12]
        Pxx /= Pxx.sum()
        entropy[i] = -np.sum(Pxx * np.log2(Pxx))
    return np.clip(entropy / np.log2(segment//2 + 1), 0, 1)

def fast_wavelet_energy(dd, level=2):
    n, m = dd.shape
    if m < 8: return np.zeros(n)
    energy = np.zeros(n)
    for i in range(n):
        try:
            wp = pywt.WaveletPacket(dd[i], 'db1', mode='symmetric', maxlevel=level)
            nodes = wp.get_level(level)
            energies = [np.sum(n.data**2) for n in nodes]
            total = sum(energies) + 1e-12
            energy[i] = energies[-1] / total
        except:
            energy[i] = 0.0
    return energy

def fast_mahal_robust(dd):
    n, m = dd.shape
    if n < 3 or m < 2: return np.zeros(n)
    try:
        if n > 50:
            idx = np.random.choice(n, 50, replace=False)
            sample = dd[idx]
        else:
            sample = dd
        lw = LedoitWolf().fit(sample)
        cov = lw.covariance_ + 1e-6 * np.eye(m)
        mean = lw.location_
        dist = np.zeros(n)
        for i in range(n):
            delta = dd[i] - mean
            dist[i] = np.sqrt(delta @ np.linalg.pinv(cov) @ delta)
        return np.clip(dist / (np.percentile(dist, 95) + 1e-12), 0, 1)
    except:
        return np.zeros(n)

def fast_kan_layer(dd):
    variances = np.var(dd, axis=1)
    return 0.1 * (variances**2 + 0.1 * np.sin(variances))

def fast_graph_degree(dd, tau=0.8):
    n = dd.shape[0]
    if n < 2: return np.zeros(n)
    norms = np.linalg.norm(dd, axis=1, keepdims=True)
    C = np.abs(dd @ dd.T) / (norms @ norms.T + 1e-12)
    np.fill_diagonal(C, 0)
    G = (C > tau).astype(int)
    deg = G.sum(axis=1)
    return 1 - deg / (deg.max() + 1e-12)

def fast_wasserstein(dd):
    n, m = dd.shape
    if n == 0: return np.zeros(0)
    template = np.mean(dd, axis=0)
    dists = np.array([wasserstein_distance(dd[i], template) for i in range(min(n, 20))])
    if n <= 20:
        return dists / (dists.max() + 1e-12)
    return np.clip(dists.mean() * np.ones(n), 0, 1)

def fast_ae_error(dd, latent=4):
    n, m = dd.shape
    if n < 2 or m < 2: return np.zeros(n)
    try:
        U, S, Vt = np.linalg.svd(dd, full_matrices=False)
        proj = U[:, :latent] @ np.diag(S[:latent]) @ Vt[:latent, :]
        err = np.sum((dd - proj)**2, axis=1)
        return np.sqrt(err) / (np.sqrt(err).max() + 1e-12)
    except:
        return np.zeros(n)

def fast_hurst(dd):
    n, m = dd.shape
    if m < 20:
        return np.full(n, 0.5)
    h = np.zeros(n)
    max_lag = min(20, m // 4)
    lags = np.arange(2, max_lag + 1)
    log_lags = np.log(lags)
    for i in range(n):
        ts = dd[i]
        tau = np.array([np.std(ts[lag:] - ts[:-lag]) for lag in lags])
        if len(tau) < 2 or np.any(tau <= 0):
            h[i] = 0.5
            continue
        log_tau = np.log(tau)
        valid = np.isfinite(log_tau)
        if valid.sum() < 2:
            h[i] = 0.5
            continue
        poly = np.polyfit(log_lags[valid], log_tau[valid], 1)
        h[i] = np.clip(poly[0], 0.1, 2.0)
    return 1 - np.abs(h - 0.5) / 0.5

def fast_kolmogorov_ensemble(dd_pre):
    if not FAST_MODE or not USE_KOLMOGOROV_ENHANCEMENTS:
        return np.zeros(dd_pre.shape[0])
    
    scores = {}
    w = KOLMOGOROV_WEIGHTS
    if w['kc'] > 0: scores['kc'] = fast_lz_complexity(dd_pre)
    if w['spectral'] > 0: scores['spectral'] = fast_spectral_entropy(dd_pre)
    if w['wavelet'] > 0: scores['wavelet'] = fast_wavelet_energy(dd_pre)
    if w['mahal'] > 0: scores['mahal'] = fast_mahal_robust(dd_pre)
    if w['kan'] > 0: scores['kan'] = fast_kan_layer(dd_pre)
    if w['graph'] > 0: scores['graph'] = fast_graph_degree(dd_pre)
    if w['wass'] > 0: scores['wass'] = fast_wasserstein(dd_pre)
    if w['ae'] > 0: scores['ae'] = fast_ae_error(dd_pre)
    if w['hurst'] > 0: scores['hurst'] = fast_hurst(dd_pre)
    
    if not scores: return np.zeros(dd_pre.shape[0])
    
    weighted = np.zeros(dd_pre.shape[0])
    total_w = sum(w.get(k, 0) for k in scores)
    for name, arr in scores.items():
        weighted += w.get(name, 0) * arr
    return weighted / (total_w + 1e-12)

# === CORE FILTERING FUNCTIONS ===
def hybrid_preprocess(dd, intensity=0.1, alpha=0.5, p_combined=None, original_shape=None):
    if not isinstance(dd, np.ndarray) or dd.ndim != 2:
        raise ValueError("Input must be a 2D NumPy array")
    n, m = dd.shape
    if n == 0 or m == 0:
        return np.ones(n), np.zeros_like(dd)

    max_abs_dd = np.max(np.abs(dd))
    dd_norm = dd / (max_abs_dd + 1e-10) if max_abs_dd > 0 else np.zeros_like(dd)

    row_means = np.mean(dd_norm, axis=1, keepdims=True)
    row_variances = np.std(dd_norm, axis=1, keepdims=True) + 1e-10
    anomaly_component = (dd_norm - row_means) / row_variances

    weights = 1.0 + alpha * (np.random.random(n) - 0.5) * intensity
    dd_preprocessed = anomaly_component * weights[:, np.newaxis]

    dd_preprocessed += np.random.normal(0, 1e-6, dd_preprocessed.shape)

    if original_shape:
        if original_shape[1] > m:
            padding = np.zeros((n, original_shape[1] - m))
            dd_preprocessed = np.hstack((dd_preprocessed, padding))
        elif original_shape[1] < m:
            dd_preprocessed = dd_preprocessed[:, :original_shape[1]]

    return weights, dd_preprocessed

def kan_univariate_layer(inputs, grid_size=5):
    """Simple Kolmogorov-Arnold univariate layer for feature mapping."""
    x = inputs.flatten()
    grid = np.linspace(x.min(), x.max(), grid_size)
    spline_coeffs = np.polyfit(x, x**2 + 0.1 * np.sin(x), 2)  # Learnable univariate (spline approximation)
    mapped = np.polyval(spline_coeffs, x).reshape(inputs.shape)
    return mapped

def ev_anomaly_detection(dd, snr_global_light, sigma2, intensity):
    if not isinstance(dd, np.ndarray) or dd.ndim != 2:
        raise ValueError("Input dd must be a 2D NumPy array")
    n, m = dd.shape
    if n == 0 or m == 0:
        return np.zeros(n), np.zeros(n), 0.0

    row_means = np.mean(np.abs(dd), axis=1)
    data_std = np.std(dd) + 1e-10

    p_anomaly = np.zeros(n)
    for i in range(n):
        row_max = np.max(np.abs(dd[i, :]) + 1e-10)
        row_signal = dd[i, :] / (row_max + 1e-10) if row_max > 0 else np.zeros_like(dd[i, :])

        row_variance = np.std(dd[i, :]).reshape(1, -1)
        kan_feature = kan_univariate_layer(row_variance).flatten()[0] * 0.1

        dot_product = np.sum(row_signal * dd[i, :]) + kan_feature
        if np.allclose(dot_product, 0, atol=1e-8):
            logging.warning(f"Zero dot_product detected for row {i}, adding perturbation.")
            dot_product += np.random.normal(0, data_std * 1e-4, 1)[0]

        norm = np.linalg.norm(dot_product) + 1e-10
        phase_angle = np.angle(dot_product / norm) if norm > 0 else 0
        phase_i = phase_angle * (3.0 + 3.0 * np.std(dd[i, :]) / (row_means[i] + 1e-10))

        anomaly_score = np.abs(phase_i) * snr_global_light / (sigma2 + 1e-10) * intensity
        if anomaly_score <= 0:
            anomaly_score = data_std * 1e-4
        p_anomaly[i] = np.clip(anomaly_score / (1 + anomaly_score), 0, 1)

    efficiency = np.mean(p_anomaly)

    anomaly_threshold = np.percentile(p_anomaly, 75)
    anomaly_map = (p_anomaly > anomaly_threshold).astype(int)

    return p_anomaly, anomaly_map, efficiency

def solve_ralf1_linear(dQ2, Q2, a, iterations=50):
    Np, Nf = dQ2.shape
    if Q2.shape[0] != Nf or a.shape[0] != Np:
        raise ValueError("Input shapes mismatch: dQ2 (Np, Nf), Q2 (Nf,), a (Np,)")

    A = np.zeros((Np, Np))
    b = np.zeros(Np)
    scale = 1e-3
    for k in range(Np):
        b[k] = np.sum(dQ2[k, :] * Q2) * scale
        for i in range(Np):
            A[k, i] = np.sum(dQ2[k, :] * dQ2[i, :]) * scale

    A += np.eye(Np) * 1e-2

    x = np.ones(Np) * 1e-6
    X_sum = np.zeros(Np)
    total_iterations = iterations
   
    for _ in range(iterations):
        indices = np.random.permutation(Np)
        x_iter = np.ones(Np) * 1e-6
        for i in indices:
            z = b[i] - np.dot(A[i, indices != i], x_iter[indices != i])
            A_ii = A[i, i] if A[i, i] != 0 else 1e-6
            x_iter[i] = z / A_ii
            x_iter[i] = np.clip(x_iter[i], -a[i], a[i])
        X_sum += x_iter
   
    x = X_sum / total_iterations
    return x

def filter_ralf1_matrix(source_matrix, intensity=1.0, iterations=50):
    Np, Nf = source_matrix.shape
    dQ2 = source_matrix
    Q2 = np.mean(source_matrix, axis=0)
    a = np.ones(Np) * np.max(np.abs(source_matrix))
    x = solve_ralf1_linear(dQ2, Q2, a, iterations)
    x_norm = x / (np.max(np.abs(x)) + 1e-10)
    x_norm = intensity * x_norm + (1 - intensity)
    filtered_matrix = source_matrix.copy()
    for i in range(Np):
        filtered_matrix[i, :] *= x_norm[i]
    return filtered_matrix

def hybrid_geo_filter(dd, intensity=0.1, iterations=15, alpha=0.5):
    if not isinstance(dd, np.ndarray) or dd.ndim != 2:
        raise ValueError("Input must be a 2D NumPy array")
    n, m = dd.shape
    if n == 0 or m == 0:
        return np.zeros_like(dd), np.zeros(n), np.zeros(n)

    dd_norm = dd / (np.max(np.abs(dd)) + 1e-10) * (10000.0 if n > 100 else 1000.0)

    weights_pre_light, dd_pre_light = hybrid_preprocess(dd_norm, intensity, alpha)

    if np.any(np.isnan(dd_pre_light)):
        logging.warning("NaN detected in dd_pre_light, replacing with zeros.")
        dd_pre_light = np.nan_to_num(dd_pre_light, nan=0.0)

    signal_norm = np.linalg.norm(dd_pre_light) + 1e-10
    noise_norm = np.linalg.norm(dd_pre_light - dd_norm) + 1e-10
    snr_global_light = 10 * np.log10(signal_norm**2 / noise_norm**2)
    sigma2 = np.var(dd_norm) + 1e-10
    p_anomaly, anomaly_map, efficiency = ev_anomaly_detection(dd_pre_light, snr_global_light, sigma2, intensity)

    if USE_KOLMOGOROV_ENHANCEMENTS and FAST_MODE:
        p_kolmo = fast_kolmogorov_ensemble(dd_pre_light)
        p_anomaly = np.maximum(p_anomaly, 0.4 * p_kolmo)

    snr_rows = np.zeros(n)
    for i in range(n):
        signal_power = np.linalg.norm(dd_pre_light[i, :])**2 + 1e-10
        noise_power = np.linalg.norm(dd_pre_light[i, :] - dd_norm[i, :])**2 + 1e-10
        snr_rows[i] = 10 * np.log10(signal_power / noise_power)
    snr_median = np.median(snr_rows)
    snr_threshold = max(np.percentile(snr_rows, 70.09), snr_median + 0.05)
    snr_anomaly_map = (snr_rows > snr_threshold).astype(int)
    p_combined = np.maximum(p_anomaly, snr_anomaly_map)

    filtered_matrix = np.zeros_like(dd_norm)
    high_prob_rows = np.where(p_combined >= 0.5)[0]
    if len(high_prob_rows) == 0:
        high_prob_rows = np.arange(n)
    if len(high_prob_rows) > 0:
        dd_high = dd_norm[high_prob_rows, :]
        weights_pre, dd_preprocessed = hybrid_preprocess(dd_high, intensity=intensity, alpha=alpha)
        logging.info(f"Rank of dd_preprocessed: {np.linalg.matrix_rank(dd_preprocessed)}/{len(high_prob_rows)}")
        filtered_matrix[high_prob_rows, :] = filter_ralf1_matrix(dd_preprocessed, intensity=intensity, iterations=iterations)
    else:
        weights_pre, dd_preprocessed = hybrid_preprocess(dd_norm[:1, :], intensity, alpha=alpha)
        filtered_matrix[0, :] = filter_ralf1_matrix(dd_preprocessed, intensity=intensity, iterations=iterations)[0, :]

    low_prob_rows = np.where(p_combined < 0.5)[0]
    if len(low_prob_rows) > 0:
        filtered_matrix[low_prob_rows, :] = dd_norm[low_prob_rows, :] * 0.1

    signal_norm = np.linalg.norm(filtered_matrix) + 1e-10
    noise_norm = np.linalg.norm(filtered_matrix - dd_norm) + 1e-10
    snr_global = 10 * np.log10(signal_norm**2 / noise_norm**2)
    p_anomaly, anomaly_map, efficiency = ev_anomaly_detection(dd_preprocessed, snr_global, sigma2, intensity)

    # Extend p_anomaly to full size
    p_anomaly_full = np.zeros(n)
    p_anomaly_full[high_prob_rows] = p_anomaly
    if len(high_prob_rows) == 0:
        p_anomaly_full[0] = p_anomaly[0]

    snr_rows = np.zeros(n)
    for i in range(n):
        signal_power = np.linalg.norm(dd_pre_light[i, :])**2 + 1e-10
        noise_power = np.linalg.norm(dd_pre_light[i, :] - dd_norm[i, :] )**2 + 1e-10
        snr_rows[i] = 10 * np.log10(signal_power / noise_power)
    snr_median = np.median(snr_rows)
    snr_threshold = max(np.percentile(snr_rows, 70.09), snr_median + 0.05)
    snr_anomaly_map = (snr_rows > snr_threshold).astype(int)
    p_combined = np.maximum(p_anomaly_full, snr_anomaly_map)

    return filtered_matrix, p_combined, snr_anomaly_map

def XFilter(dd0, key3=1):
    # stddQ2=np.std(dd0)
    # srdQ2=np.mean(dd0,axis=0)
    # dd0=dd0-srdQ2
    if not key3==1:
       dd0=XFilterB(dd0)        
    dd, p_combined, snr_anomaly_map = hybrid_geo_filter(dd0, intensity=0.1, iterations=15, alpha=0.5)   
    # dd=(dd-np.mean(dd,axis=0))*stddQ2/np.std(dd)+srdQ2
    return dd  # Return the filtered matrix

# === Faddeev Residual Refinement ===
def faddeev_residual_refinement(dQ4_B, dQ4_A):
    return hybrid_geo_filter(dQ4_B - dQ4_A, 0.15, 10, 0.6)[0]



# =============================================================================
# END OF CODE
# =============================================================================

import psutil
def CheckTemp(aTemp=MaxTemp):
    battery = psutil.sensors_battery()
    if not battery is None:
        return int(battery.power_plugged)
    else:
        return 1

DETERM=0.1

def RandomQ(Nfx,NQRandm_=0):
    while not CheckTemp():
        tm.sleep(60)
    QRandm_=np.asarray(range(512),float)
    NQRandm=0  
    KK=3e6
    liiX=np.zeros(Nfx,float)
    pp=0
    while pp<0.7:
        for ii in range(3):
            try:                
                z=(QRandm_[NQRandm]+1)/KK           
                atim0=tm.time()        
                tm.sleep(z) 
                atim=tm.time()     
                dd=int(((atim-atim0)/z-1)/1000)
                zz=np.asarray(sample(list(range(Nfx)),Nfx),float)/KK
                LFib1340(dd).shuffle(zz)  
                LFib1340(int(2*dd/(1+np.sqrt(5)))).shuffle(QRandm_)
                
                if NQRandm>0:
                    liiX=liiX+zz
                NQRandm=NQRandm+1
            except:
                NQRandm=0                

        k2, pp = scp.skewtest(liiX)
            
    rr2=[[],[]]
    rr2[0]= liiX.copy()
    rr2[1]= np.asarray(np.linspace(0,Nfx-1,Nfx),int)
    m=[[rr2[j][l] for j in range(len(rr2))] for l in range(len(rr2[0]))]         
    m.sort(key=itemgetter(0))                  
    rr2=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]  
    liiXX=np.asarray(rr2[1],int)
    return liiXX

def filterFourierQ(arxx,arb,NNew,NChan,key=-1,key2=0): 
    if key2==-2:
        arxx=np.concatenate((arxx, arxx[::-1]))
        arb=np.concatenate((arb, arb[::-1]))
        Nnl=2*NNew  
    else:
        Nnl=NNew
    Nfl=int(len(arb)/NChan)
    Nfl_=int(len(arxx)/NChan)
    
    #ar_=np.zeros(Nnl,float)
    farx=np.zeros(8*Nnl,float)-1e32
    
    az=int(np.floor(Nfl/Nnl))
    
    gg0=0    
    for l in range(NChan):        
        for i in range(az):
            ar_=arb[Nfl-(az-i)*Nnl+Nfl*l:Nfl+Nnl-(az-i)*Nnl+Nfl*l].copy()
            gg0=gg0+np.sum(ar_*ar_)
            #ar_=ar_[::-1].copy()
            ar_=ar_-ar_[0]
            #ar_=ar_-np.mean(ar_)
            ar_x=ar_[::-1].copy()
            ar_=np.concatenate((ar_x,-2*((key<0)-.5)*ar_))
            ar_=np.concatenate((ar_[::-1],ar_))
            ar_=ar_-ar_[len(ar_)-1]
            ar_=np.concatenate((ar_,-ar_))
            ar__=abs(np.fft.fft(ar_))
            farx=np.maximum(farx,ar__)
    #gg0=np.sqrt(gg0)/(NChan*az*Nnl)
     
    gg=0
    arxr=arxx.copy()
    for l in range(NChan):      
        ar_=arxx[Nfl_*(l+1)-Nnl:Nfl_*(l+1)].copy()
        #par_= savgol_filter(ar_, 14, 5)
        par_=ar_.copy()
        #ar_=ar_-ar_[0]
        ar_=ar_[::-1].copy()
        ar_=ar_-par_[0]
        #ar_=ar_-ar_[len(ar_)-1]
        #ar_=ar_-np.mean(ar_)
        ar_x=ar_[::-1].copy()
        ar_=np.concatenate((ar_,-2*((key<0)-.5)*ar_x))
        ar_=np.concatenate((ar_[::-1],ar_))
        ar_=ar_-ar_[len(ar_)-1]
        ar_=np.concatenate((ar_,-ar_))
        farxx=np.fft.fft(ar_)    
        mfarxx=np.abs(farxx)+1e-32  
        farxxx=np.zeros(8*Nnl,complex)    

        A1=np.exp(DETERM*np.std(np.log(mfarxx[5:8*Nnl-4]))+np.mean(np.log(mfarxx[5:8*Nnl-4])))
 
        for j in range(1,8*Nnl):
            if mfarxx[j]>A1:
                farxxx[j]=(farxx[j]/mfarxx[j])*np.sqrt(farx[j]*mfarxx[j]) 
                #farxxx[j]=(farxx[j]/mfarxx[j])*farx[j]
                       
        farxxx[0]=0*farxx[0]          
                   
        aaa=np.fft.ifft(farxxx).real
        aaa=aaa[3*Nnl:5*Nnl].copy()
        aaa=-((aaa[0:Nnl]-2*((key<0)-.5)*aaa[2*Nnl:Nnl-1:-1])/2)#aaa[0:Nnl][::-1]#
        arxr[Nfl_-Nnl+Nfl_*l:Nfl_+Nfl_*l]=aaa.copy()
        #arxr[Nfl_-Nnl+Nfl_*l:Nfl_+Nfl_*l]= savgol_filter(arxr[Nfl_-Nnl+Nfl_*l:Nfl_+Nfl_*l], 14, 5)
        gg=gg-arxr[Nfl_-Nnl+Nfl_*l]+arb[Nfl_-Nnl+Nfl_*l-1]

    gg=gg/NChan  
    if key2==-2:
        aln=int(len(arxr)/2)
        arxr=arxr[:aln]
    return arxr

def RALF1Calculation2(arr_bx,arr_c,Nf,NNew,NNew0,NChan,Nhh,iProc,Nproc):
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    sz=Nf*NChan
    MM=3
    Nzz=2#Nhh
    Nhh=Nzz
    
    Ndel=MM
    NCh=int(np.ceil(sz/Ndel))  
    Ndel0=MM
    NCh0=int(np.ceil(sz/Ndel0))   
          
    arr_b=np.asarray(arr_bx,float)
    
    hh=0  
    hh_=0      
    Nzz0=Nzz
    Atim_0=tm.time()  
    dd1_x=[]
    dd2_x=[]  
    while hh<Nhh:  
        if hh==0:  
            AMX=np.zeros((Nhh,sz),float)
            AMN=np.zeros((Nhh,sz),float)
            dd10=np.zeros((Nhh+1,sz),float)
            dd20=np.zeros((Nhh+1,sz),float)
            max_dd1=np.zeros((Nhh+1,sz),float)-1e32
            min_dd2=np.zeros((Nhh+1,sz),float)+1e32            
            arr_bbbxxx1=np.zeros((Nhh,sz),float)
            arr_bbbxxx2=np.zeros((Nhh,sz),float)

            Nzz=Nzz0
            rr2=np.zeros((Nhh+1,sz),float)
            rr2[hh]=np.asarray(arr_b,float)
        if hh==1:
            Nzz=int(Nzz0/2+1)
            
        wwwww=0
        while not wwwww:
            try:
                arrrxx=hkl.load("ralfrez.rlf2")
                wwwww=1
            except:
                tm.sleep(0.1)

        if len(arrrxx)>=Nproc: 
            return rr2[hh]/0 
        
        liix=np.zeros((sz,sz),int) 
        dQ3=np.zeros((sz,sz),float)
        mDD=np.zeros((sz,sz),float)  
        
        aa=RandomQ(sz)  
        liiB=np.concatenate((aa,aa,aa))  
        aa=RandomQ(sz) 
        liiD=np.concatenate((aa,aa,aa))          

        rrr=rr2[hh].copy()
        # rrr_=rrr[0]
        # rrr[1:]=np.diff(rrr)
            
        for i in range(sz):    
            liix[i]=liiB[liiD[i]:sz+liiD[i]].copy()#range(sz)#
            dQ3[i]=rrr[liix[i]].copy()

        astart=np.inf
        dQ3=dQ3.reshape((sz*sz))
        astart=dQ3[0]            
        dQ3[1:]=np.diff(dQ3)
        dQ3[0]=0
        dQ3=dQ3.reshape((sz,sz))
        dQ3_=dQ3.copy()
        
        D=np.std(dQ3)
        
        R4=np.ones(sz,float)  
        for l in range(NChan):
            R4[Nf-NNew+Nf*l:Nf+Nf*l]=0
            
        for i in range(sz):     
            mDD[i]=R4[liix[i]].copy() 
        
        seqq=(dQ3.reshape(sz*sz))[1:]*np.ceil(0.5*(1/(mDD.reshape(sz*sz)==1)[0:sz*sz-1]+1/(mDD.reshape(sz*sz)==1)[1:]))
        seqq_=  mDD.reshape(sz*sz)*0
        seqq_[0]=1
        seqq_[1:]=(abs(seqq)==np.inf)+np.isnan(seqq)
        seqq_=seqq_.reshape((sz,sz))
        seqq=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqq)),float) 
        seqq=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqq)),float)

        ##########################################       
        aa=RandomQ(sz) 
        NumFri0=np.concatenate((aa, aa, aa))  
        aa=RandomQ(sz) 
        NumFri0_=np.concatenate((aa, aa, aa)) 
        aa=RandomQ(sz) 
        rR0=np.concatenate((aa, aa, aa))   
        aa=RandomQ(sz) 
        liiC=np.concatenate((aa, aa, aa)) 
        aa=RandomQ(sz)
        aa=aa-np.mean(aa)
        r5=aa.copy()
        r5=20*r5*D/np.std(r5)
        r5=r5-np.mean(r5)
        r5=np.concatenate((r5, r5))
        aa=RandomQ(sz) 
        ss4=np.concatenate((aa, aa, aa, aa))                         
        zz=0  
        WW=0   
                
        while zz<Nzz and WW>-2*Nhh: 
            KY=-1
            # dQ3=dQ3_.copy()   
            dQ3mx=np.zeros((sz,sz))-1e32
            dQ3mn=np.zeros((sz,sz))+1e32
            dQ3num=np.zeros((sz,sz),int)
            NumFri=NumFri0_[NumFri0[ss4[zz]]:NumFri0[ss4[zz]]+2*sz].copy()
            NumFri_=NumFri0[NumFri0_[ss4[zz]]:NumFri0_[ss4[zz]]+2*sz].copy()
            while np.min(dQ3num)<2:                  
                aa=ss4[0]
                ss4[0:len(ss4)-1]=ss4[1:].copy()
                ss4[len(ss4)-1]=aa
                rR=rR0[liiC[ss4[zz]]:liiC[ss4[zz]]+2*sz].copy()
                rR_=rR0[liiC[ss4[len(ss4)-2*sz+zz]]:liiC[ss4[len(ss4)-2*sz+zz]]+2*sz].copy()
                
                kk=-1
                xxx=0
                while kk <Ndel-1 and xxx==0:  
                    kk=kk+1                      
                    ii=int(kk*NCh)
                    k=-1
                    while k<Ndel0-1 and xxx==0:     
                        k=k+1                       
                        i=int(k*NCh0) 
                        dQ4=np.zeros((NCh,NCh0),float)
                        mDD4=np.zeros((NCh,NCh0),float)
                        mDD4_=np.zeros((NCh,NCh0),float)
                        mDD4_A=np.zeros((NCh,NCh0),float) 
                        mDD4_B=np.zeros((NCh,NCh0),float)                                    
                        for ll in range(NCh0):
                            dQ4[:,ll]=dQ3[NumFri[ii:ii+NCh],NumFri_[i+ll]].copy()
                            mDD4[:,ll]=mDD[NumFri[ii:ii+NCh],NumFri_[i+ll]].copy()
       
                        dQ4=dQ4-np.mean(dQ4)
                        P=np.zeros(3,float)
                        
                        for ll in range(NCh0):
                            vvv=r5[rR[ss4[ll]+zz]:rR[ss4[ll]+zz]+NCh]*(1-mDD4[:,ll])
                            mDD4_A[:,ll]=vvv.copy()
                            mDD4_B[:,ll]=vvv.copy()#vvv[::-1].copy()
                       
                        mDD4_B=mDD4_A.copy()#*0
                        # KY=KY*(-1)
                        # mDD4_A=KY*mDD4_A*(mDD4_A>0)
                        # mDD4_B=KY*mDD4_B*(mDD4_B>0)
                        
                        # nNxA=sum(sum(mDD4==1))     
                        # nNxA_=sum(sum(mDD4==0))  
                        # if nNxA>nNxA_ and nNxA_>0:  
                        # if nNxA_>0:  
                        seqA0=(dQ4.reshape(NCh*NCh0))[1:]*np.ceil(0.5*(1/(mDD4.reshape(NCh*NCh0)==1)[0:NCh*NCh0-1]+1/(mDD4.reshape(NCh*NCh0)==1)[1:]))
                        seqA0_=  mDD4.reshape(NCh*NCh0)*0
                        seqA0_[0]=1
                        seqA0_[1:]=(abs(seqA0)==np.inf)+np.isnan(seqA0)
                        seqA0_=seqA0_.reshape((NCh,NCh0))
                        seqA0=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqA0)),float) 
                        seqA0=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqA0)),float)

                        
                        seqA=(dQ4.reshape(NCh*NCh0))[1:]*np.ceil(0.5*np.fabs(1/(1*(mDD4.reshape(NCh*NCh0)==1)[0:NCh*NCh0-1]-1*(mDD4.reshape(NCh*NCh0)==1)[1:])))
                        seqA_=  mDD4.reshape(NCh*NCh0)*0
                        seqA_[0]=1
                        seqA_[1:]=(abs(seqA)==np.inf)+np.isnan(seqA)
                        seqA_=seqA_.reshape((NCh,NCh0))
                        seqA=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqA)),float) 
                        seqA=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqA)),float)

                        # dQ4_B=   XFilterB(mDD4_A+dQ4)-((mDD4_A))#,0)
                        dQ4_A= XFilterB( mDD4_A+dQ4)#((mDD4_A)) #,0)#,1)#
                        # eeB=XFilterB(eeB)
                        dQ4_B= ( mDD4_A+dQ4)#-((mDD4_A)) #,0)#,1)#  
                        # dQ4_B=XFilterB(dQ4_B)
                        # dQ4_A=  -XFilterB(mDD4_B-dQ4)+((mDD4_B))#,0)
                        # dQ4_A=XFilterB(dQ4_A)
                        # dQ4_B=(dQ4_B+dQ4_A)/2
                        # # asr1=abs(dQ4_B-xxxx)>abs(dQ4_A-xxxx)
                        # # asr2=abs(dQ4_B-xxxx)<abs(dQ4_A-xxxx)                    
                        # # dQ4_B=dQ4_B*asr1+dQ4_A*asr2+(dQ4_B+dQ4_A)*(asr1==asr2)/2    

                        # # dQ4_B = faddeev_residual_refinement(dQ4_B, dQ4_A)                           
                        # # # # # dQ4_B=dQ4_B-dQ4_A
                        # dQ4_A=dQ4_B.copy()
                        
                        seqB=((dQ4_A).reshape(NCh*NCh0))[1:]*np.ceil(0.5*np.fabs(1/(1*(mDD4.reshape(NCh*NCh0)==1)[0:NCh*NCh0-1]-1*(mDD4.reshape(NCh*NCh0)==1)[1:])))
                        seqB=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqB)),float) 
                        seqB=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqB)),float)
                        seqC=((dQ4_B).reshape(NCh*NCh0))[1:]*np.ceil(0.5*np.fabs(1/(1*(mDD4.reshape(NCh*NCh0)==1)[0:NCh*NCh0-1]-1*(mDD4.reshape(NCh*NCh0)==1)[1:])))
                        seqC=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqC)),float) 
                        seqC=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqC)),float)
                        try:    
                            P_1=P.copy()
                            P_2=P.copy()
                            P_1[0:2]=np.polyfit(seqA,seqB,1)
                            P_2[0:2]=np.polyfit(seqA,seqC,1)
                            # print(100*scp.pearsonr(seqA,seqB)[0])
                            # P_1[0]=np.std(seqB)/np.std(seqA)
                            # P_1[1]=np.mean(seqB)-P_1[0]*np.mean(seqA)  
                            # P_1[2]=0
                            # P_2[0]=np.std(seqC)/np.std(seqA)
                            # P_2[1]=np.mean(seqC)-P_2[0]*np.mean(seqA)  
                            # P_2[2]=0
                            if P_1[0]>0 and P_2[0]>0:# and 100*scp.pearsonr(seqA,seqB)[0]>0 and 100*scp.pearsonr(seqA,seqC)[0]>0:# and not (abs(P_1[0]-1)>0.5 or abs(P_2[0]-1)>0.5):   
                                dQ4_A=(dQ4_A-P_1[1])/P_1[0]
                                dQ4_B=(dQ4_B-P_2[1])/P_2[0]                                     
                                # dQ4_B=faddeev_residual_refinement(dQ4_B,dQ4_A)
                                dQ4_B=XFilter(dQ4_B-dQ4_A)                                
                                dQ4_A=dQ4_B.copy() 
                                
                                seqD=((dQ4_B).reshape(NCh*NCh0))[1:]*np.ceil(0.5*np.fabs(1/(1*(mDD4.reshape(NCh*NCh0)==1)[0:NCh*NCh0-1]-1*(mDD4.reshape(NCh*NCh0)==1)[1:])))
                                seqD=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqD)),float) 
                                seqD=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqD)),float)
                                # # P_1[0:2]=np.polyfit(seqA,seqB,1)
                                # # dQ4_B=(dQ4_B-P_1[1])/P_1[0]
                                # # dQ4_A=dQ4_B.copy() 
                                # P_1[0]=np.std(seqD)/np.std(seqA)
                                # P_1[1]=np.mean(seqD)-P_1[0]*np.mean(seqA)  
                                # P_1[2]=0
                                # P_2[0]=np.std(seqD)/np.std(seqA)
                                # P_2[1]=np.mean(seqD)-P_2[0]*np.mean(seqA)  
                                # P_2[2]=0
                                # dQ4_A=(dQ4_A-P_1[1])/P_1[0]
                                # dQ4_B=(dQ4_B-P_2[1])/P_2[0]  
                                zcoef=100*scp.pearsonr(seqA,seqD)[0]
                                # print(zcoef)
                                # if zcoef>0:
                                #     P_1[0:2]=np.polyfit(seqA,seqD,1)
                                #     P_2[0:2]=np.polyfit(seqA,seqD,1)
                                #     if P_1[0]<0 and P_2[0]<0:
                                #         dQ4_A=(dQ4_A-P_1[1])/P_1[0]
                                #         dQ4_B=(dQ4_B-P_2[1])/P_2[0] 
                                for ll in range(NCh0):                                        
                                    seqA=dQ3_[NumFri[ii:ii+NCh],NumFri_[i+ll]]+seqA0_[:,ll]*np.maximum(dQ3mx[NumFri[ii:ii+NCh],NumFri_[i+ll]]-dQ3_[NumFri[ii:ii+NCh],NumFri_[i+ll]],dQ4_A[:,ll]-dQ3_[NumFri[ii:ii+NCh],NumFri_[i+ll]])
                                    seqB=dQ3_[NumFri[ii:ii+NCh],NumFri_[i+ll]]+seqA0_[:,ll]*np.minimum(dQ3mn[NumFri[ii:ii+NCh],NumFri_[i+ll]]-dQ3_[NumFri[ii:ii+NCh],NumFri_[i+ll]],dQ4_B[:,ll]-dQ3_[NumFri[ii:ii+NCh],NumFri_[i+ll]])
                                    ab=dQ3num[NumFri[ii:ii+NCh],NumFri_[i+ll]]*dQ3mx[NumFri[ii:ii+NCh],NumFri_[i+ll]]
                                    ba=dQ3num[NumFri[ii:ii+NCh],NumFri_[i+ll]]*dQ3mn[NumFri[ii:ii+NCh],NumFri_[i+ll]]
                                    dQ3mx[NumFri[ii:ii+NCh],NumFri_[i+ll]]=(ab+seqA)/(dQ3num[NumFri[ii:ii+NCh],NumFri_[i+ll]]+1)
                                    dQ3mn[NumFri[ii:ii+NCh],NumFri_[i+ll]]=(ba+seqB)/(dQ3num[NumFri[ii:ii+NCh],NumFri_[i+ll]]+1)
                                    
                                    dQ3num[NumFri[ii:ii+NCh],NumFri_[i+ll]]=dQ3num[NumFri[ii:ii+NCh],NumFri_[i+ll]]+1
                                    dQ3[NumFri[ii:ii+NCh],NumFri_[i+ll]]=(dQ3mx[NumFri[ii:ii+NCh],NumFri_[i+ll]]+dQ3mn[NumFri[ii:ii+NCh],NumFri_[i+ll]])/2

                                ll=ll+1
                            # else:
                            #     xxx=1
                            
                        except:
                            xxx=1
                        # else:     
                        #     xxx=1

            # dQ3=dQ3_.copy()

            if np.isnan(np.std(dQ3mx)) or np.isnan(np.std(dQ3mn)):
                xxx=1
            if xxx==0:  
                # dQ3mxL=dQ3mx.copy()
                # dQ3mnL=dQ3mn.copy()                 
                if zz==0:
                    AsrXMx=dQ3mx.copy()
                    AsrXMn=dQ3mn.copy()      
                else:  
                    # AsrXMx=(AsrXMx*zz+dQ3mx)/(zz+1)
                    # AsrXMn=(AsrXMn*zz+dQ3mn)/(zz+1)
                    AsrXMx=np.maximum(AsrXMx,dQ3mx)
                    AsrXMn=np.minimum(AsrXMn,dQ3mn)
                    # AsrXMx=(AsrXMx*zz+np.maximum(AsrXMx,dQ3mx))/(zz+1)
                    # AsrXMn=(AsrXMn*zz+np.minimum(AsrXMn,dQ3mn))/(zz+1)

                dQ3=(AsrXMx+AsrXMn)/2
                AsrXMx_=AsrXMx.copy()
                AsrXMn_=AsrXMn.copy()
                
                WW=0                                    
                zz=zz+1
            else:
                aa=ss4[0]
                ss4[0:len(ss4)-1]=ss4[1:].copy()
                ss4[len(ss4)-1]=aa
                WW=WW-1                      
                
        wwwww=0
        while not wwwww:
            try:
                arrrxx=hkl.load("ralfrez.rlf2")
                wwwww=1
            except:
                tm.sleep(0.1)

        if len(arrrxx)>=Nproc: 
            return rr2[hh]/0 
        
        hh0=hh
        if not WW<0:   
            1==1
            #AsrXMx_=filter_ralf1_matrix(AsrXMx_, intensity=1.0, iterations=100)
            seqq_1=((AsrXMx_).reshape(sz*sz))[1:]*np.ceil(0.5*(1/(mDD.reshape(sz*sz)==1)[0:sz*sz-1]+1/(mDD.reshape(sz*sz)==1)[1:]))
            seqq_1=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqq_1)),float) 
            seqq_1=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqq_1)),float)
            
            #AsrXMn_=filter_ralf1_matrix(AsrXMn_, intensity=1.0, iterations=100)
            seqq_2=((AsrXMn_).reshape(sz*sz))[1:]*np.ceil(0.5*(1/(mDD.reshape(sz*sz)==1)[0:sz*sz-1]+1/(mDD.reshape(sz*sz)==1)[1:]))
            seqq_2=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqq_2)),float) 
            seqq_2=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqq_2)),float)
            P_1[0:2]=np.polyfit(seqq,seqq_1,1)
            P_2[0:2]=np.polyfit(seqq,seqq_2,1)
            if 100*scp.pearsonr(seqq,seqq_1)[0]>0 and 100*scp.pearsonr(seqq,seqq_2)[0]>0 and P_1[0]>0 and P_2[0]>0:#not abs(P[0]-1)>0.5 and :
                1==1
                AsrXMx_=(AsrXMx_-P_1[1])/P_1[0] 
                AsrXMn_=(AsrXMn_-P_2[1])/P_2[0] 
            else: 
                WW=-1
        if not WW<0:         
            # dd0=0
            # asr1=abs(AsrXMx_-dd0)>abs(AsrXMn_-dd0)
            # asr2=abs(AsrXMx_-dd0)<abs(AsrXMn_-dd0)                    
            # AsrXMx_=AsrXMx_*asr1+AsrXMn_*asr2+(AsrXMx_+AsrXMn_)*(asr1==asr2)/2
            # AsrXMn_=AsrXMx_.copy()

            if not astart==np.inf:
                dd1=AsrXMx_.reshape((sz*sz))
                dd1=astart+np.cumsum(dd1)
                AsrXMx_=dd1.reshape((sz,sz))
                dd1=AsrXMn_.reshape((sz*sz))
                dd1=astart+np.cumsum(dd1)
                AsrXMn_=dd1.reshape((sz,sz))
                  
            AsrXMx_=(AsrXMx_+AsrXMn_)/2 
            AsrXMn_=AsrXMx_.copy()
            
            aMx_=0
            aMn_=0
            aMx=np.zeros(sz,float)-1e32
            aMn=np.zeros(sz,float)+1e32
            for i in  range(sz):
                aMx[liix[i]]=np.maximum(aMx[liix[i]],AsrXMx_[i])
                aMn[liix[i]]=np.minimum(aMn[liix[i]],AsrXMn_[i])
                # aMx_=(aMx_*i+aMx)/(i+1)#
                # aMn_=(aMn_*i+aMn)/(i+1)#
            aMx_=aMx.copy()                    
            aMn_=aMn.copy()
            ann=sum(np.isnan(aMx_ + aMn_))
            if ann==0:               
                if hh==0: 
                    AMX[hh]=aMx_.copy()
                    AMN[hh]=aMn_.copy()  
                    arr_bbbxxx1[hh]=AMX[hh].copy()
                    arr_bbbxxx2[hh]=AMN[hh].copy()
                    D1=10
                else:
                    AMX[hh]=aMx_.copy()
                    AMN[hh]=aMn_.copy()
                    arr_bbbxxx1[hh]=(arr_bbbxxx1[hh-1]*(hh-1)+AMX[hh])/(hh)
                    arr_bbbxxx2[hh]=(arr_bbbxxx2[hh-1]*(hh-1)+AMN[hh])/(hh)  
                    dd1_=[]
                    dd2_=[]                    
                    D1_=[]
                    D2_=[]
                    for l in range(NChan):
                        dd1_.append(arr_bbbxxx1[hh-1][Nf-NNew+Nf*l:Nf+Nf*l].copy())  
                        dd2_.append(arr_bbbxxx2[hh-1][Nf-NNew+Nf*l:Nf+Nf*l].copy())  
                        D1_.append((arr_bbbxxx1[hh]-arr_bbbxxx1[hh-1])[Nf-NNew+Nf*l:Nf+Nf*l].copy())
                        D2_.append((arr_bbbxxx2[hh]-arr_bbbxxx2[hh-1])[Nf-NNew+Nf*l:Nf+Nf*l].copy())
                    D1_=np.asarray(D1_,float)
                    D2_=np.asarray(D2_,float)
                    dd1_=np.asarray(dd1_,float)
                    dd2_=np.asarray(dd2_,float)                  
                    D1=np.std(D1_+D2_)/np.std(dd1_+dd2_)*np.sqrt(hh+1)/np.sqrt(1+(hh+1))  
                    
                ann=1  
                try:
                    # for l in range(NChan):                             
                    #     rrr[Nf*l:Nf+Nf*l]= savgol_filter(rrr[Nf*l:Nf+Nf*l], 14, 5)

                    dd1=AMX[hh].copy()#(filterFourierQ(AMX[hh],rrr,NNew,NChan))
                    dd2=AMN[hh].copy()#(filterFourierQ(AMN[hh],rrr,NNew,NChan)) 
                    # dd0=(dd1+dd2)/2
                    # asr1=abs(dd1-dd0)>abs(dd2-dd0)
                    # asr2=abs(dd1-dd0)<abs(dd2-dd0)                    
                    # dd1=dd1*asr1+dd2*asr2+(dd1+dd2)*(asr1==asr2)/2
                    # dd2=dd1.copy()
                    if sum(np.abs(dd1+dd2)==np.inf)==0 and D1>DETERM:  
                        sr2_1=[]
                        sr2_2=[]
                        sarr_c=[]
                        for l in range(NChan):  
                            sr2_1.append((dd1)[Nf*l:Nf-NNew0+Nf*l])
                            sr2_2.append((dd2)[Nf*l:Nf-NNew0+Nf*l])
                            sarr_c.append(rrr[Nf*l:Nf-NNew0+Nf*l])
                        sr2_1=np.asarray(sr2_1,float)
                        sr2_2=np.asarray(sr2_2,float)
                        sarr_c=np.asarray(sarr_c,float) 
                        sr2_1=sr2_1.reshape((len(sr2_1)*len(sr2_1[0])))
                        sr2_2=sr2_2.reshape((len(sr2_2)*len(sr2_2[0])))
                        sarr_c=sarr_c.reshape((len(sarr_c)*len(sarr_c[0])))   
                    
                        P_1=P.copy()
                        P_2=P.copy()
                        
                        # sr2_1=(sr2_1+sr2_2)/2 
                        # sr2_2=sr2_1.copy()
                        P_1[0:2]=np.polyfit(sarr_c,sr2_1,1)
                        P_2[0:2]=np.polyfit(sarr_c,sr2_2,1) 
                        P_1[0]=np.std(sr2_1)/np.std(sarr_c)
                        P_1[1]=np.mean(sr2_1)-P_1[0]*np.mean(sarr_c) 
                        P_2[0]=np.std(sr2_2)/np.std(sarr_c)
                        P_2[1]=np.mean(sr2_2)-P_2[0]*np.mean(sarr_c) 
                        
                        #abs(P_1[0]-P_2[0])/(abs(P_1[0])+abs(P_2[0]))<0.2
                        if abs(P_1[0]-P_2[0])/(abs(P_1[0])+abs(P_2[0]))<1 and P_1[0]>0 and P_2[0]>0:# and 100*scp.pearsonr(sarr_c,((sr2_1-P_1[1])/P_1[0]+(sr2_2-P_2[1])/P_2[0]))[0]>0:                         
                            for l in range(NChan):  
                                dd1=(dd1-P_1[1])/P_1[0]
                                dd2=(dd2-P_2[1])/P_2[0]
                            max_dd1[hh]=rr2[hh].copy()
                            min_dd2[hh]=rr2[hh].copy()
                            for l in range(NChan):
                                if hh==0:
                                     max_dd1[hh,Nf-NNew+Nf*l:Nf+Nf*l]=dd1[Nf-NNew+Nf*l:Nf+Nf*l].copy()
                                     min_dd2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=dd2[Nf-NNew+Nf*l:Nf+Nf*l].copy()
                                else:
                                      max_dd1[hh,Nf-NNew+Nf*l:Nf+Nf*l]=np.maximum(max_dd1[hh-1,Nf-NNew+Nf*l:Nf+Nf*l],
                                                                        (dd1[Nf-NNew+Nf*l:Nf+Nf*l]))
                                      min_dd2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=np.minimum(min_dd2[hh-1,Nf-NNew+Nf*l:Nf+Nf*l],
                                                                        (dd2[Nf-NNew+Nf*l:Nf+Nf*l]))
                                                                               
                            hh=hh+1
                            ann=0         
                            # rr2[hh]=(max_dd1[hh-1]+min_dd2[hh-1])/2
                            
                            dd1=max_dd1[hh-1].copy()
                            dd2=min_dd2[hh-1].copy()
                            dd0=(dd1+dd2)/2
                            asr1=abs(dd1-dd0)>abs(dd2-dd0)
                            asr2=abs(dd1-dd0)<abs(dd2-dd0)                    
                            rr2[hh]=dd1*asr1+dd2*asr2+(dd1+dd2)*(asr1==asr2)/2
                            rr2[hh]=(rr2[hh-1]*(hh-1)+filterFourierQ(rr2[hh],rrr,NNew,NChan))/hh 
                            
                            ##rr2[hh]=(rr2[hh-1]*(hh-1)+rr2[hh])/hh                             
                            sr2=[]
                            sarr_c=[]
                            for l in range(NChan):  
                                sr2.append(rr2[hh,Nf-NNew+int(0.3*(NNew-NNew0))+Nf*l:Nf-NNew0+Nf*l].copy())
                                sarr_c.append(arr_c[(NNew-NNew0)*l+int(0.3*(NNew-NNew0)):NNew-NNew0+(NNew-NNew0)*l].copy())
                            sr2=np.asarray(sr2,float)
                            sarr_c=np.asarray(sarr_c,float) 
                            sr2=sr2.reshape((len(sr2)*len(sr2[0])))
                            sarr_c=sarr_c.reshape((len(sarr_c)*len(sarr_c[0])))   
                            
                            P[0:2]=np.polyfit(sarr_c,sr2,1)                              
                            
                            if P[0]>0:# and 100*scp.pearsonr(sarr_c,sr2)[0]>0:#and :
                                P[0]=np.std(sr2)/np.std(sarr_c)
                                P[1]=np.mean(sr2)-P[0]*np.mean(sarr_c) 
                                for l in range(NChan):  
                                    rr2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=(rr2[hh,Nf-NNew+Nf*l:Nf+Nf*l]-P[1])/P[0]
                                    rr2[hh,Nf*l:Nf-NNew+Nf*l]=arr_b[Nf*l:Nf-NNew+Nf*l].copy()
                                    # rr2[hh,:]= savgol_filter(rr2[hh,:], 14, 5)
                                #P[0:2]=np.polyfit(sarr_c,sr2,1)   
                                
                                #D1=0.2
                                if hh==Nhh:
                                    if 1:#abs(P[0]-1.)<0.5 and sum(abs(rr2[hh])==np.inf)==0 and D1<1: 
                                        anamef="fralf.tmp"
                                        fo = open(anamef, "w")
                                        Atim_1=tm.time()   
                                        fo.write(str(iProc)+' it=%d, K=%3.2f Tm=%5.2f\n'%(hh,D1,Atim_1-Atim_0))
                                        fo.close()
                                        return np.asarray(rr2[hh],float)
                                    else:
                                        hh=hh0
                            else:                                
                                hh=hh0
                        else:
                            hh=hh0
                except:
                    hh=hh0

        if hh0==hh:            
            hh_=hh_+1
            if hh>1:
                hh=hh-2
            else:
                hh=0
                
            if hh_>2*Nhh:
                return rr2[hh]/0                                                                                             

def RALf1FiltrQ(*args):  
    pid = win32api.GetCurrentProcessId()
    handle = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS, True, pid)
    win32process.SetPriorityClass(handle, priorityclasses[1])
    
    ########################### 
    www=0
    while (not www):
        try:
            with open("filtrx.tmp", "rb") as file:
                args=file.read()
                file.close() 
                www=1
        except:
            tm.sleep(0.1)
    [NChan, Nf, NNew, Nhh, Nproc]=np.asarray(struct.unpack('5i', args[0:20]),np.int32)
    anamef="maxmin_.xxx"
    fo = open(anamef, "w")
    fo.write('%s '%(Nproc)+'\n')
    fo.close() 
    arr_bb=np.asarray(struct.unpack('%sf'%(Nf), args[20:(4*Nf+20)]),np.float32)
    ###########################
      
    NNew0=int(NNew*1.2) 
    
    #Nhh=1
    Nhh0=Nhh
    astar0=np.inf
    # astar0=arr_bb[0]

    # arr_bb[0]=0
    
    arr_b=np.asarray(arr_bb,float)
    
    arr_c=[]
    aa=RandomQ(NNew0)
    aa=aa-np.mean(aa)
    #aa= savgol_filter(aa, 14, 5)
    aa=20*aa/np.std(aa)*np.std(arr_b)
    ss4=np.concatenate((aa, aa, aa, aa)) #*0    
    for l in range(NChan):
        arr_b[Nf-NNew0+Nf*l:Nf+Nf*l]=(astar0==np.inf)*arr_b[Nf-NNew0+Nf*l-1]+ss4[l:NNew0+l]
        arr_b[Nf*l:Nf+Nf*l]=np.cumsum(arr_b[Nf*l:Nf+Nf*l])
        arr_bb[Nf*l:Nf+Nf*l]=np.cumsum(arr_bb[Nf*l:Nf+Nf*l])
        arr_c.append(arr_b[Nf-NNew0+Nf*l:Nf-NNew+Nf*l].copy()) 
        
    arr_c=np.asarray(arr_c,float)
    arr_c=arr_c.reshape((len(arr_c)*len(arr_c[0])))
    
    while 1==1: 
        hh=0
        ann=0
        arr_bbx=np.zeros((Nhh*2,Nf),float)
        Koef=np.zeros(Nhh0+1,float)                
        KoefA=np.zeros(Nhh0+1,float)
        NumIt=3*(int((Nhh0)/2)+1)

        while hh<Nhh:
            if hh<Nhh:    
                arr_bbbxxx=RALF1Calculation2(arr_b,arr_c,Nf,NNew0,NNew,NChan,NumIt,args[0],Nproc)                
                
                wwwww=0
                while not wwwww:
                    try:
                        arrrxx=hkl.load("ralfrez.rlf2")
                        wwwww=1
                    except:
                        tm.sleep(0.1)

                if len(arrrxx)>=Nproc: 
                    return arr_bbbxxx/0 
                if (sum(np.abs(arr_bbbxxx)==np.inf)==0 and sum(np.isnan(arr_bbbxxx))==0):
                    Nf_=NNew0+int(NNew0*1.2)
                    if Nf_>=Nf:
                        Nf_=Nf-1
                    arr_bbbxxx_=np.zeros(Nf_*NChan,float)
                    arr_c_=[]
                    aa=RandomQ(NNew0)
                    aa=aa-np.mean(aa)
                    aa=20*aa/np.std(aa)*np.std(arr_bbbxxx)
                    ss4=np.concatenate((aa, aa, aa, aa)) 
                    for l in range(NChan):
                        arr_bbbxxx_[Nf_*l:Nf_+Nf_*l]=np.asarray(arr_bbbxxx[Nf*(l+1):Nf-Nf_-1+Nf*l:-1],float)
                        arr_c_.append(arr_bbbxxx_[Nf_-NNew0+Nf_*l:Nf_-NNew+Nf_*l].copy())
                        arr_bbbxxx_[Nf_-NNew0+Nf_*l:Nf_+Nf_*l]=(astar0==np.inf)*arr_bbbxxx_[Nf_-NNew0+Nf_*l-1]+ss4[l:NNew0+l]                     
                    arr_c_=np.asarray(arr_c_,float)
                    arr_c_=arr_c_.reshape((len(arr_c_)*len(arr_c_[0]))) 
                    arr_bbbxxx_y=RALF1Calculation2(arr_bbbxxx_,arr_c_,Nf_,NNew0,NNew,NChan,NumIt,args[0],Nproc)
                    if (sum(np.abs(arr_bbbxxx_y)==np.inf)==0 and sum(np.isnan(arr_bbbxxx_y))==0): 
                        for l in range(NChan):
                            if not astar0==np.inf:
                                if l==0:
                                    mm1=np.cumsum(arr_bb[Nf-NNew0-(Nf_-NNew0)+Nf*l:Nf*(l+1)-(Nf_-NNew0)])
                                    mm2=np.cumsum(arr_bbbxxx_y[Nf_*(l+1):Nf_-NNew0-1+Nf_*l:-1])
                                else:
                                    mm1=np.concatenate((mm1,np.cumsum(arr_bb[Nf-NNew0-(Nf_-NNew0)+Nf*l:Nf*(l+1)-(Nf_-NNew0)])))
                                    mm2=np.concatenate((mm2,np.cumsum(arr_bbbxxx_y[Nf_*(l+1):Nf_-NNew0-1+Nf_*l:-1])))                
                            else:
                                if l==0:
                                    mm1=arr_bb[Nf-NNew0-(Nf_-NNew0)+Nf*l:Nf*(l+1)-(Nf_-NNew0)].copy()
                                    mm2=arr_bbbxxx_y[Nf_*(l+1):Nf_-NNew0-1+Nf_*l:-1].copy()
                                else:
                                    mm1=np.concatenate((mm1,arr_bb[Nf-NNew0-(Nf_-NNew0)+Nf*l:Nf*(l+1)-(Nf_-NNew0)].copy()))
                                    mm2=np.concatenate((mm2,arr_bbbxxx_y[Nf_*(l+1):Nf_-NNew0-1+Nf_*l:-1].copy()))                

                        ann=(sum(np.abs(mm1)==np.inf)>0 + sum(np.isnan(mm1))>0+
                              sum(np.abs(mm2)==np.inf)>0 + sum(np.isnan(mm2))>0)
                        
                        if ann==0 and len(mm1)>1 and len(mm1)==len(mm2):                             
                            if np.std(mm1)>0 and np.std(mm2)>0:
                                coef=100*(scp.pearsonr(mm1,mm2)[0])
                                if coef>0:
                                    anamef="fralf_.tmp"
                                    fo = open(anamef, "w")
                                    fo.write(str(args[0])+' %s'%(coef)+'\n')
                                    fo.close() 
                                
                                KoefA[hh]=coef
         
                                #mm1=mm1*np.std(mm2)/np.std(mm1)                       
                                Koef[hh]=-np.std(mm1-mm2)
                                if not astar0==np.inf:
                                    arr_bbx[hh]=astar0+np.cumsum(arr_bbbxxx)
                                else:
                                    arr_bbx[hh]=arr_bbbxxx.copy()
                                #hh=hh+1
                
            # if hh==Nhh:
                                
                                try:
                                    wwwww=0
                                    while not wwwww:
                                        try:
                                            arrrxx=hkl.load("ralfrez.rlf2")
                                            wwwww=1
                                        except:
                                            tm.sleep(0.1)

                                    if len(arrrxx)>=Nproc: 
                                        return np.asarray(arr_bbx[hh],float)
                                    else:
                                        print('%s'%(coef))
                                        if coef>0:   
                                            zresu=arr_bbx[hh].copy()
                                            zresu[1:]=np.diff(zresu)
                                            zresu[0]=0
                                            arrrxx.append(np.asarray(zresu,float))
                                            wwwww=0
                                            while not wwwww:
                                                try:
                                                    hkl.dump(arrrxx,"ralfrez.rlf2")
                                                    wwwww=1
                                                except:
                                                    tm.sleep(0.6)
                                except:
                                    arrrxx=[]

#https://query1.finance.yahoo.com/v7/finance/download/LRC-USD?period1=1635554377&period2=1667097577&interval=1d&events=history&includeAdjustedClose=true

api_key = 'ONKTYPV6TAMZK464' 
 
interv="15min"
interv="Daily"

#INTRADAY
#d_intervals = {"1min","5min","15min","30min","60min"}

Lengt0=250
Ngroup=3
Nproc=Ngroup*(os.cpu_count())
Lo=1  
lSrez=0.99
aTmStop=1
NIt=3
NIter=100
DT=0.3
dNIt=4
aDecm=1
andecim=1
KPP=0

aKEY=0
def decimat(adat_):
    if Lo:
        adat_ = np.log(adat_)
    k = 0
    adat__ = np.zeros(int(len(adat_) / aDecm), float)
    for i in range(int(len(adat_) / aDecm)):
        adat__[k] = np.median(adat_[i * aDecm:i * aDecm + aDecm])
        k = k + 1
        
    arrray=adat__[1:len(adat__)].copy()
    
    indices = np.arange(len(arrray))
    non_nan_mask = ~np.isnan(arrray)
    nan_mask = np.isnan(arrray)
    
    interpolated_values = np.interp(indices[nan_mask], indices[non_nan_mask], arrray[non_nan_mask])
    
    arrray[nan_mask] = interpolated_values
    
    if Lo:
        return np.exp(arrray)
    else:
        return (arrray)
    
def fig2img ( fig ):
    fig.savefig(wrkdir +ticker+'dynamic.png',dpi=150,transparent=False,bbox_inches = 'tight')
    frame=Image.open(wrkdir +ticker+'dynamic.png')
    # fig.canvas.draw()
    # frame=Image.frombytes('RGB', fig.canvas.get_width_height(),
    #                        fig.canvas.tostring_rgb())
    return frame

def loaddata(aLengt, ticker1, key):
    adat_ = []
    with open(wrkdir + ticker1 + '.csv', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        dat = []
        i = 0
        for row in spamreader:
            if i > 0:
                dat.append(row[0])
            i = i + 1
    dat = dat[len(dat)::-1]
    try:
        dat = np.asarray(dat, float)
        arrr = np.asarray(dat, float)
    except:
        dat = dat
        arrr = np.zeros(2, int)
        
    if andecim>=1:
        sz=len(arrr)
        xold=np.asarray(range(sz),int)
        f=interp.interp1d(xold, arrr,'cubic')
        xnew=np.asarray(range((sz-1)*andecim),int)/andecim
        arrr=f(xnew)
        sz=len(arrr)
        arrr=arrr[sz-aLengt*aDecm:sz]
        # %varexp --plot arrr 
        
    return arrr, adat_

try:
    nnams_=hkl.load(wrkdir +"name.rlf1")
    ii=len(nnams_)
except:
    WhO = [
        'XRP', "AAVE","ANKR","AVAX","AXS","BAL","BTC","BCH","ADA","CGLD","LINK","COMP",
        "CRV","ICP","MANA","DOGE","EOS","ETH","ETC","GRT","LTC","LRC","MKR","NU",
        "OMG","ORN","MATIC","SHIB","SKL","SOL","XLM","SNX","XTZ","UNI"
    ]

    WhO= ['MANA', 'ANKR', 'XLM', 'ORN', 'COMP', 'AVAX', 'NU']
    
    def getcsv(WhO,xYears,wrkdir):
        for i in range(len(WhO)):
            nm=WhO[i]
            anm = "https://www.cryptodatadownload.com/cdd/Binance_%sUSDT_d.csv"%(WhO[i])
            try:
                df = pd.read_csv(anm, delimiter=",", skiprows=[0]) 
                # Extract 'Close' column as a float-type NumPy array
                arrr = df['Close']
                arrr.to_csv(wrkdir+nm+'.csv', index=False)
            except:
                anm=anm
                
    if GetSCV:        
        getcsv(WhO,5,wrkdir)

    arrrxxR = []
    nams = []
    lenar = []

    for i in range(len(WhO)):
        i0 = 0
        while i0 < 6 and not i0 < 0:
            try:
                arrrxx1, adat1_ = loaddata(Lengt0, WhO[i], 1)
                arrrxxR.append(arrrxx1)
                lenar.append(len(arrrxx1))
                nams.append(WhO[i])
                i0 = -1
            except:
                i0 = i0 + 1
                #tm.sleep(10)

    llar = int(np.median(np.asarray(lenar, int)))
    nnams_ = []
    aaer = []
    for i in range(len(nams)):
        if len(arrrxxR[i]) >= llar:
            aer = decimat(arrrxxR[i])
            aaer.append(aer[len(aer) - int(llar / aDecm) + 1:].copy())
    aaer = np.asarray(aaer, float)
            
    arrrxxR_ = []
    ii = 0
    for i in range(len(nams)):
        if len(arrrxxR[i]) >= llar:
            arrrxxR_.append(aaer[ii])
            nnams_.append(nams[ii])
            ii = ii + 1
    arrrxxR_ = np.asarray(arrrxxR_, float)
               
    aKEY=1
    for uuii in range(len(nnams_)):               
        aname=nnams_[uuii]
        ticker=aname+"YLLL"
        ticker1=aname
        ticker2=aname
        try:                
            hkl.dump(arrrxxR_[uuii],wrkdir + aname+"dat.rlf1")
        except:
            os.mkdir(wrkdir)
            hkl.dump(arrrxxR_[uuii],wrkdir + aname+"dat.rlf1")
    hkl.dump(nnams_,wrkdir + "name.rlf1")

def RALF1Cella(*arrgs_): 
    pid = win32api.GetCurrentProcessId()
    handle = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS, True, pid)
    win32process.SetPriorityClass(handle, priorityclasses[1])
    
    www=0
    while (not www):
        try:
            with open("maxmindat.tmp", "rb") as file:
                bin_arrgsx=file.read()
                file.close()     
            arrgs=bin_arrgsx            
            [Nf,anI,NNew,hhhx,aMM,aNN,Numproc]=np.asarray(struct.unpack('7i', arrgs[0:28]),np.int32)
            # anamef="maxmin_.tmp"
            # fo = open(anamef, "w")
            # fo.write('%s '%(Numproc)+'\n')
            # fo.close() 
            binZDat=np.asarray(struct.unpack('%sf'%(anI*Nf), arrgs[28:(4*anI*Nf+28)]),np.float32)
            ZDat=binZDat.reshape((anI,Nf))
            ss0=np.asarray(struct.unpack('%si'%(int((len(arrgs)-(4*anI*Nf+28))/4)), arrgs[(4*anI*Nf+28):]),np.int32)

            www=1
        except:
            tm.sleep(0.1)
    ww=1      
    #ZData_=np.mean(ZDat[:,0:Nf-NNew],axis=0)
    # aa=RandomQ(Nf)                        
    # ss0=np.concatenate((aa, aa, aa))
    while (ww):
        dd=ZDat.copy()
        for ii in range(len(dd)):
            dd[ii]=np.cumsum(dd[ii])  
        mdd4=dd*0                       
        aa=RandomQ(Nf)                        
        ss4=np.concatenate((aa, aa, aa))
        liix=np.zeros((anI,Nf),int)
        mdd4_=np.zeros(Nf,float)
        mdd4_[0:Nf-NNew]=1
        for i in range(anI):  
            liix[i]=ss4[ss0[i+hhhx+ww:i+Nf+hhhx+ww]].copy()
            dd[i]=(dd[i])[liix[i]].copy() 
            mdd4[i]=mdd4_[liix[i]].copy()
        
        astart=np.inf
        dd=dd.reshape((anI*Nf))  
        astart=dd[0]
        dd[1:]=np.diff(dd)
        dd[0]=0                        
        dd=dd.reshape((anI,Nf))                                 
        D=np.std(dd)
        dd0=dd.copy()    
                                
        P=np.zeros(3,float)
        PP=1
        
        seq0=(dd.reshape(len(dd)*len(dd[0])))[1:]*np.ceil(0.5*(1/(mdd4.reshape(len(dd)*len(dd[0]))==1)[0:len(dd)*len(dd[0])-1]+1/(mdd4.reshape(len(dd)*len(dd[0]))==1)[1:]))
        seq0=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seq0)),float) 
        seq0=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seq0)),float)    
        
        dd_AA=dd*0
        dd_BB=dd*0
        dd_CC=dd*0
        dd_Num=dd*0
        dd_CCmx=dd*0-1e32
        dd_CCmn=dd*0+1e32


        while np.min(dd_Num)<2:    
            aa=RandomQ(Nf)  
            aa=aa-np.mean(aa)                 
            ss4_=np.concatenate((aa, aa, aa))                                      
            DD_=[]
            for hhhc in range(anI):
                vvv=(ss4_[hhhc:hhhc+Nf]*(1-mdd4[hhhc])).copy()
                DD_.append(vvv[::-1].copy())
            DD_=np.asarray(DD_,float)                              
            DD_=(DD_/np.std(DD_))*D*20

            #DD_=DD_*0
            KY=-1
            for ii in range(aNN):   
                for jj in range(aMM):    
                    dd1=dd[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)].copy()
                    dd1=dd1-np.mean(dd1)
                    mdd4_=mdd4[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)].copy()
                    
                    seqA0=(dd1.reshape(len(dd1)*len(dd1[0])))[1:]*np.ceil(0.5*(1/(mdd4_.reshape(len(dd1)*len(dd1[0]))==1)[0:len(dd1)*len(dd1[0])-1]+1/(mdd4_.reshape(len(dd1)*len(dd1[0]))==1)[1:]))
                    seqA0_=  mdd4_.reshape(len(dd1)*len(dd1[0]))*0
                    seqA0_[0]=1
                    seqA0_[1:]=(abs(seqA0)==np.inf)+np.isnan(seqA0)
                    seqA0_=seqA0_.reshape((len(dd1),len(dd1[0])))
                    seqA0=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqA0)),float) 
                    seqA0=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqA0)),float)    
                    
                    seqA=(dd1.reshape(len(dd1)*len(dd1[0])))[1:]*np.ceil(0.5*np.fabs(1/(1*(mdd4_.reshape(len(dd1)*len(dd1[0]))==1)[0:len(dd1)*len(dd1[0])-1]-1*(mdd4_.reshape(len(dd1)*len(dd1[0]))==1)[1:])))                                                                                                       
                    seqA_=  mdd4_.reshape(len(dd1)*len(dd1[0]))*0
                    seqA_[0]=1
                    seqA_[1:]=(abs(seqA)==np.inf)+np.isnan(seqA)
                    seqA_=seqA_.reshape((len(dd1),len(dd1[0])))
                    seqA=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqA)),float) 
                    seqA=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqA)),float)    
                    
                    DD__A=DD_[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)].copy()
                    DD__B=DD__A.copy()#
                    # KY=KY*(-1)
                    # DD__A=KY*DD__A*(DD__A>0)
                    # DD__B=KY*DD__B*(DD__B>0)
                    
                    if len(dd1)>1 and len(dd1[0])>=len(dd1):
                        eeA= XFilterB(DD__A+dd1)#Filter( DD__A+dd1)#-((DD__A)) #,0)#,1)#
                        # eeB=XFilterB(eeB)
                        eeB= ( DD__A+dd1)#-((DD__A)) #,0)#,1)#    
                        # eeA=XFilterB(eeA)
                        # eeB=(eeB+eeA)/2
                        # # asr1=abs(eeB-xxxx)>abs(eeA-xxxx)
                        # # asr2=abs(eeB-xxxx)<abs(eeA-xxxx)                    
                        # # eeB=eeB*asr1+eeA*asr2+(eeB+eeA)*(asr1==asr2)/2    
                        
                        # # # # eeB = faddeev_residual_refinement(eeB, eeA)                           
                        # # # # # eeB=eeB-eeA
                        # eeA=eeB.copy()
                        
                        dd_AA[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]=(eeA).copy()#*(eeB>0)*((eeA+eeB)>0)
                        dd_BB[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]=(eeB).copy()#*(eeA<0)*((eeA+eeB)<0)
                    dd2_1=(dd_AA)[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)].copy()
                    seqB=(dd2_1.reshape(len(dd2_1)*len(dd2_1[0])))[1:]*np.ceil(0.5*np.fabs(1/(1*(mdd4_.reshape(len(dd2_1)*len(dd2_1[0]))==1)[0:len(dd2_1)*len(dd2_1[0])-1]-1*(mdd4_.reshape(len(dd2_1)*len(dd2_1[0]))==1)[1:])))
                    seqB=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqB)),float) 
                    seqB=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqB)),float)    
                    dd2_2=(dd_BB)[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)].copy()
                    seqC=(dd2_2.reshape(len(dd2_2)*len(dd2_2[0])))[1:]*np.ceil(0.5*np.fabs(1/(1*(mdd4_.reshape(len(dd2_2)*len(dd2_2[0]))==1)[0:len(dd2_2)*len(dd2_2[0])-1]-1*(mdd4_.reshape(len(dd2_2)*len(dd2_2[0]))==1)[1:])))
                    seqC=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqC)),float) 
                    seqC=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqC)),float)    
    
                    try:
                        P_1=P.copy()
                        P_2=P.copy()
                        P_1[0:2]=np.polyfit(seqA,seqB,1)
                        P_2[0:2]=np.polyfit(seqA,seqC,1)

                        # print(100*scp.pearsonr(seqA,seqB)[0])
                        # P_1[0]=np.std(seqB)/np.std(seqA)
                        # P_1[1]=np.mean(seqB)-P_1[0]*np.mean(seqA)  
                        # P_1[2]=0
                        # P_2[0]=np.std(seqC)/np.std(seqA)
                        # P_2[1]=np.mean(seqC)-P_2[0]*np.mean(seqA)  
                        # P_2[2]=0
                        if P_1[0]>0 and P_2[0]>0:# and 100*scp.pearsonr(seqA,seqB)[0]>0 and 100*scp.pearsonr(seqA,seqC)[0]>0:
                            eeA=(eeA-P_1[1])/P_1[0]
                            eeB=(eeB-P_2[1])/P_2[0]   
                            eeA=XFilter(eeB-eeA)
                            eeB=eeA.copy()
                                                        
                            seqD=(eeB.reshape(len(dd2_1)*len(dd2_1[0])))[1:]*np.ceil(0.5*np.fabs(1/(1*(mdd4_.reshape(len(dd2_1)*len(dd2_1[0]))==1)[0:len(dd2_1)*len(dd2_1[0])-1]-1*(mdd4_.reshape(len(dd2_1)*len(dd2_1[0]))==1)[1:])))
                            seqD=np.asarray(list(filter(lambda x: abs(x)!= np.inf, seqD)),float) 
                            seqD=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqD)),float)    
                            zcoef=100*scp.pearsonr(seqA,seqD)[0]
                            
                            # P_1[0]=np.std(seqD)/np.std(seqA)
                            # P_1[1]=np.mean(seqD)-P_1[0]*np.mean(seqA)  
                            # P_1[2]=0
                            # P_2[0]=np.std(seqD)/np.std(seqA)
                            # P_2[1]=np.mean(seqD)-P_2[0]*np.mean(seqA)  
                            # P_2[2]=0
                            # eeA=(eeA-P_1[1])/P_1[0]
                            # eeB=(eeB-P_2[1])/P_2[0]  
                            
                            # print(zcoef)
                            # if zcoef>0:
                            aanum=dd_Num[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)].copy()
                            seqAA=dd0[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]+np.maximum(dd_CCmx[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]-dd0[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)],(dd2_1)-dd0[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)])
                            seqBB=dd0[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]+np.minimum(dd_CCmn[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]-dd0[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)],(dd2_2)-dd0[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)])
                            dd_CCmx[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]=(dd_CCmx[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]*aanum+np.maximum((dd_CCmx[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]),seqAA))/(aanum+1)
                            dd_CCmn[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]=(dd_CCmn[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]*aanum+np.minimum((dd_CCmn[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]),seqBB))/(aanum+1)                               
                            dd[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]=(dd_CCmx[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]+dd_CCmn[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)])/2
                            dd_Num[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]=aanum+1
                    except:
                        PP=0
        # dd_CCmx=filter_ralf1_matrix(dd_CCmx, intensity=1.0, iterations=100)
        # dd_CCmn=filter_ralf1_matrix(dd_CCmn, intensity=1.0, iterations=100)                
        if not PP==0:               
            dd0=(dd_CCmx+dd_CCmn)/2
            asr1=abs(dd_CCmx-dd0)>abs(dd_CCmn-dd0)
            asr2=abs(dd_CCmx-dd0)<abs(dd_CCmn-dd0)                    
            dd_CCmx=dd_CCmx*asr1+dd_CCmn*asr2+(dd_CCmx+dd_CCmn)*(asr1==asr2)/2
            dd_CCmn=dd_CCmx.copy()                                                            
            if not astart==np.inf: 
                # dd_AA=dd_CC*PP*(dd_CC>0)
                # dd_BB=dd_CC*PP*(dd_CC<0)     
                dd_AA=dd_CCmx.reshape((anI*Nf))
                dd_AA=astart+np.cumsum(dd_AA)
                dd_AA=dd_AA.reshape((anI,Nf))
                dd_BB=dd_CCmn.reshape((anI*Nf))
                dd_BB=astart+np.cumsum(dd_BB)
                dd_BB=dd_BB.reshape((anI,Nf))
                
            dd_AA=(dd_AA+dd_BB)
            dd_BB=dd_AA.copy()  
            aMx=np.zeros(Nf,float)-1e32
            aMn=np.zeros(Nf,float)+1e32
            for i in range(anI):
                aMx[liix[i]]=np.maximum(aMx[liix[i]],dd_AA[i])
                aMn[liix[i]]=np.minimum(aMn[liix[i]],dd_BB[i])

            aRes=np.asarray([aMx, aMn],np.float32)     
            
            aRes[0][1:]=np.diff(aRes[0])
            aRes[0][0]=0              
            aRes[1][1:]=np.diff(aRes[1])  
            aRes[1][0]=0                     
            ww=0
            
        www=0
        while (not www):
            try:
                with open("maxmin.tmp", "rb") as file:                    
                    areza=file.read()
                    file.close()   
                www=1

            except:
                tm.sleep(0.1)
            
        nsz=int(len(areza)/2/4/Nf)
        if nsz<Numproc:
            www=0
            while (not www):
                try:   
                    areza=np.asarray(struct.unpack('%sf'%(int(len(areza)/4)), areza),np.float32)
                    areza=np.asarray(areza.reshape((nsz,2,Nf)),np.float32)
                    areza_=np.zeros((nsz+1,2,Nf),np.float32)
                    areza_[0:nsz]=areza.copy()
                    areza_[nsz]=aRes.copy()
                        
                    with open("maxmin.tmp", "wb") as file:
                        ssssss=np.asarray(areza_.reshape(((nsz+1)*2*Nf)),np.float32)     
                        sssss=struct.pack("%sf"%(len(ssssss)),*ssssss)                                          
                        file.write(sssss)                        
                        file.close() 
                    www=1
                except:
                    tm.sleep(0.1)  
            ww=ww+1  
        else:
            ww=0
            www=0
            while (not www):
                try:
                    # anamef="maxmin_x.tmp"
                    # fo = open(anamef, "w")
                    # fo.write('%s '%(aRes)+'\n')
                    # fo.close() 
                    www=1
                except:
                    tm.sleep(0.1)  
                    
            break
            #os._exit(os.EX_OK)
    return aRes



warnings.filterwarnings("ignore", category=RuntimeWarning) 
#nnams__=nnams_.copy()
nnams_=hkl.load(wrkdir + "name.rlf1")
Nii=len(nnams_)
if __name__ == '__main__': 
    for uuii in range(Nii):
        nnams_=hkl.load(wrkdir + "name.rlf1")
        Nii=len(nnams_)
        aname=nnams_[uuii]
        try:    
            os.remove(wrkdir + aname+".ralf")
        except:
            1==1
        print(aname)
        ticker=aname+"YLLL"
        ticker1=aname
        ticker2=aname
        try:
            dill.load_session(wrkdir + aname+".ralf")
            nnams_=hkl.load(wrkdir + "name.rlf1")
            Nii=len(nnams_)
        except:
            ImApp=[]
            try:
                arrrxx=hkl.load(wrkdir + aname+"dat.rlf1")
            except:
    
                if not ticker1==ticker2:
                    arrrxx1,adat1_=loaddata(Lengt0,ticker1,aKEY)
                    arrrxx1=np.asarray(arrrxx1,float)
                    arrrxx2,adat2_=loaddata(Lengt0,ticker2,aKEY)
                    arrrxx2=np.asarray(arrrxx2,float)
                    lnm=min(len(arrrxx1),len(arrrxx2))
                    arrrxx=arrrxx1[len(arrrxx1)-lnm:]/arrrxx2[len(arrrxx2)-lnm:]
                else:
                    ticker=ticker1
                    arrrxx,adat_=loaddata(Lengt0,ticker,aKEY)
                    arrrxx=np.asarray(arrrxx,float)
                
                arrrxx=decimat(arrrxx)
                    
                try:                
                    hkl.dump(arrrxx,wrkdir + aname+"dat.rlf1")
                except:
                    os.mkdir(wrkdir)
                    hkl.dump(arrrxx,wrkdir + aname+"dat.rlf1")
                
            esz=len(arrrxx)
            arr_rezDzRez=[[] for j in range(esz)]
            kkk=0
            out=0   
    
            arrr=np.asarray(arrrxx).copy()  
    
            arrr=np.asarray(arrr,float)    
            Lengt=len(arrr)
            Nf=Lengt
            
            nn=int(Nf*DT)             
            NNew=int(Nf*0.5)  
            Nf=Nf+nn        
            ar0=np.asarray(arrr[0:])           
            
            arr_z=np.zeros(Nf,float)
            arr_z[0:Nf-NNew]=arrr[0:Nf-NNew].copy()
            arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]
              
            adat0=''
                                   
            Arr_AAA=np.zeros((NIter*Nproc,Nf),float) 
            arr_rezBz=np.zeros(Nf,float)
            arr_rezBzz=arr_rezBz.copy()
            
            all_rezAz=np.zeros((NIter,Nf),float)
            arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]  
            all_RezN=np.zeros((Ngroup,NIter,Nf),float)+np.inf
            all_RezM=np.zeros((Ngroup,NIter,Nf),float)-np.inf
            dd1a=np.zeros((Ngroup,NIter,Nf),float)-np.inf
            dd2a=np.zeros((Ngroup,NIter,Nf),float)+np.inf
            all_RezNM=np.zeros((Ngroup,NIter,Nf),float)
            all_RezMM=np.zeros((Ngroup,NIter,Nf),float)
            all_RezNM_=np.zeros((Ngroup,NIter,Nf),float)
            all_RezMM_=np.zeros((Ngroup,NIter,Nf),float)

            #argss=[[0] for j in range(Nproc)]    
    
            hh0=0
            hhh=0
            hhh_=0
                    
            Koef_=[]
            ZZ=0
            key=0
            try:
                dill.load_session(wrkdir + aname+".ralf")
                nnams_=hkl.load(wrkdir + "name.rlf1")
                Nii=len(nnams_)
            except:    
                fig = plt.figure()
                axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
                axes_ = fig.add_axes([0, 0, 0.3, 0.3])     
                axes__ = fig.add_axes([0.4, 0, 0.3, 0.3])
                
                if Lo:
                    axes.semilogy(ar0, 'r.')
                    axes.semilogy(arr_z, 'go-')  #cut data used for model
                    axes.grid(True, which="both", ls="-")
                else:
                    axes.plot(ar0, 'r.')
                    axes.plot(arr_z, 'go-')  #cut data used for model
                    axes.grid(True, which="both", ls="-")
    
                axes.text(4, 4, 'Reflection on Action of Lorentz Forces-1, #2011612714  \n\n Course = %s, start = %s, step = %s * %s / %s'%(aname,adat0,interv,aDecm,andecim),
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes_.transAxes,color='blue', fontsize=14)        
    
                try:
                    frame=fig2img(fig)  
                except:
                    frame=fig2img(fig) 
                ImApp.append(frame)
                cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)        
                gray_sz1=len(cimg[0])
                gray_sz2=len(cimg)
                aDur=2
                fourcc = cv.VideoWriter_fourcc(*'MP4V')
                out = cv.VideoWriter(wrkdir + aname+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
                for icl in range(len(ImApp)):
                    cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
                    out.write(cimgx[0:gray_sz2,0:gray_sz1,:]) 
                out.release()
                plt.show()
                del(out)
                
        #key=13
        while hhh_<aTmStop and not key == 13: 
            Aprocess=[]
            if hhh==int(NIter/5):
                if hhh_<aTmStop-1:
                    try:
                        os.remove(wrkdir + aname+".rlf1")
                    except:
                        hhh_=hhh_
                    nnn=int(nn*0.5)
                    aTmStop=6
                    hh0=0
                    hhh=0
                    arr_z=np.zeros(Nf+nnn,float)
                    arr_z[0:nnn]=arr_rezBz[0:nnn].copy()
                    arr_z[nnn:Nf]=arr_rezBz[nnn:Nf].copy()
                    arr_z[Nf:]=arr_z[Nf-1]
                    Lengt=Lengt+nnn
                    ar0=arr_z[0:Lengt-nnn].copy()
                    Nf=Nf+nnn
                    arr_rezBz=np.zeros(Nf,float)
                    arr_rezBzz=arr_rezBz.copy()
                    arr_rezMx=  np.zeros((Ngroup,Nf),float)
                    arr_rezMn=  np.zeros((Ngroup,Nf),float)
                    Arr_AAA=np.zeros((NIter*Nproc,Nf),float) 
                    all_rezAz=np.zeros((NIter,Nf),float)
                    arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]  
                    all_RezN=np.zeros((Ngroup,NIter,Nf),float)+1e32
                    all_RezM=np.zeros((Ngroup,NIter,Nf),float)-1e32
                    dd1a=np.zeros((Ngroup,NIter,Nf),float)-1e32
                    dd2a=np.zeros((Ngroup,NIter,Nf),float)+1e32
                    all_RezNM=np.zeros((Ngroup,NIter,Nf),float)
                    all_RezMM=np.zeros((Ngroup,NIter,Nf),float)
                    all_RezNM_=np.zeros((Ngroup,NIter,Nf),float)
                    all_RezMM_=np.zeros((Ngroup,NIter,Nf),float)
                    hhh_=hhh_+1
                else:
                    hhh_=hhh_+1
                    ZZ=1
            
            if ZZ==0:                  
                try:
                    [hhha,Arr_AAA]=(hkl.load(wrkdir + aname+".rlf1"))       
                except:            
                    hhha=hh0-1
                               
                if hh0>=hhha: 
                    if Lo:
                        arr_A=np.log(arr_z) 
                    else:
                        arr_A=arr_z.copy()
                        
                    Asr=np.mean(arr_A[0:Nf-NNew])
                    arr_A=arr_A-Asr
                    Klg=np.power(10,np.floor(np.log10(np.max(abs(arr_A)))))
                    arr_A=arr_A/Klg
                    # x=np.asarray(range(len(arr_A)),float)
                    # P=np.polyfit(x[:Nf-NNew],arr_A[:Nf-NNew],1)   
                    # arr_A=arr_A-np.polyval(P, x)
                    # arr_A[Nf-NNew:]=arr_A[Nf-NNew-1]
                    
                    program =wrkdir + "RALF1FiltrX_lg.py"
                    NChan=1
                    ##########################
                    ssssss=np.asarray([NChan, Nf, NNew, NIt, Nproc],np.int32)                    
                    sssss=np.asarray(arr_A,np.float32)                        
                    argss=(struct.pack("%si"%(len(ssssss)),*ssssss)+
                                struct.pack("%sf"%(Nf),*sssss))
                    with open("filtrx.tmp", "wb") as file:
                        file.write(argss)
                        file.close()  
                    ##########################
                    
                    # for iProc in range(Nproc):
                    #     argss[iProc]=["%s"%iProc, "%s"%NChan, "%s"%NNew, "%s"%NIt]#"%s"%(iProc+1)]
                    #     for i in range(Nf):
                    #         argss[iProc].append(str("%1.6f"%(arr_A[i])))   
                    #     argss[iProc].append("%d"%Nproc)                    
                    
                    try:
                        wwwww=1
                        while wwwww and wwwww<100:
                            try:
                                arezAMx_=hkl.load("ralfrez.rlf2")
                                wwwww=0
                            except:
                                tm.sleep(0.1)
                                wwwww=wwwww+1
                                
                        if wwwww>0:
                            arezAMx_=[]

                        if len(arezAMx_)==Nproc or arezAMx_==[]:
                            wwwww=1
                            while wwwww and wwwww<100:
                                try:
                                    hkl.dump([],"ralfrez.rlf2")
                                    wwwww=0
                                except:
                                    tm.sleep(0.6)
                                    wwwww=wwwww+1

                    except:
                            wwwww=0
                            while not wwwww:
                                try:
                                    hkl.dump([],"ralfrez.rlf2")
                                    wwwww=1
                                except:
                                    tm.sleep(0.6)
                    
                    arezAMx_=[] 
                    # for iProc in range(Nproc):
                    #     aaa=RALf1FiltrQ(zip(argss))
                    #     arezAMx_.append(aaa)
                    
                    pool = mp.Pool(processes=Nproc)
                    try:
                        #pool.map(RALf1FiltrQ, argss)
                        pool.map_async(RALf1FiltrQ,zip(argss))
                        wwwwww=0
                        while not wwwwww:
                            wwwww=0
                            while not wwwww:
                                try:
                                    arezAMx_=hkl.load("ralfrez.rlf2")
                                    wwwww=1
                                except:
                                    tm.sleep(0.1)
                            nsz=len(arezAMx_)
                            if nsz>=Nproc:
                                wwwwww=1 
                                pool.terminate()
                            else:
                                tm.sleep(0.1)
                    except:
                        arezAMx_=hkl.load("ralfrez.rlf2")
                    #arezAMx= np.asarray(arezAMx,float)[0,:,:]
                    del(pool)
                    
                    if arezAMx_==[]:                        
                        wwwww=0
                        while not wwwww:
                            try:
                                arezAMx_=hkl.load("ralfrez.rlf2")
                                wwwww=1
                            except:
                                tm.sleep(0.1)
                                
                    if len(arezAMx_)>0:                        
                        wwwww=0
                        while not wwwww:
                            try:
                                hkl.dump(arezAMx_,"ralfrez_.rlf2")
                                wwwww=1
                            except:
                                tm.sleep(0.6)
                        
                    arezAMx= np.asarray(arezAMx_,float)                    
                    
                    wwwww=0
                    while not wwwww:
                        try:
                            hkl.dump([],"ralfrez.rlf2")
                            wwwww=1
                        except:
                            tm.sleep(0.6)
                    # arezAMx= (np.asarray(arezAMx,float)+np.polyval(P, x))*Klg+Asr
                    
                    arezAMx= (np.asarray(arezAMx,float))*Klg+Asr
                     
                    for iGr in range(Ngroup):
                        Arr_AAA[iGr*NIter*int(Nproc/Ngroup)+hh0*int(Nproc/Ngroup):iGr*NIter*int(Nproc/Ngroup)+(hh0+1)*int(Nproc/Ngroup)]=(
                            arezAMx[int(iGr*(Nproc/Ngroup)):int((iGr+1)*(Nproc/Ngroup))]).copy()
                                                 
                    hkl.dump([hh0+1,Arr_AAA], wrkdir + aname+".rlf1")  
                    [hhha,Arr_AAA]=(hkl.load(wrkdir + aname+".rlf1"))
                
                WrtTodr=1
                if hhh>=hhha-1:   
                    #WrtTodr=1
                    aDur=4
                                    
                aNN=3
                aMM=3

                for iGr in range(Ngroup):  
                    ZDat=(Arr_AAA[iGr*NIter*int(Nproc/Ngroup)+((hh0))*int(Nproc/Ngroup):iGr*NIter*int(Nproc/Ngroup)+(hh0+min(dNIt, hhha-hh0))*int(Nproc/Ngroup)]).copy()
                    
                    # ZDat=Arr_AAA[iGr*NIter*int(Nproc/Ngroup)+max(0,0)*int(Nproc/Ngroup):iGr*NIter*int(Nproc/Ngroup)+((hhha))*int(Nproc/Ngroup)]ZDat=np.exp(Arr_AAA[iGr*NIter*int(Nproc/Ngroup)+max(0,(hh0+1)-dNIt)*int(Nproc/Ngroup):iGr*NIter*int(Nproc/Ngroup)+(hh0+1)*int(Nproc/Ngroup)])
                
                    if iGr==0:
                        xxxx=ZDat.copy()
                    else:
                        xxxx=np.concatenate((xxxx, ZDat))
                ZDat=xxxx.copy()   
                anI=len(ZDat)
                # for i in range(anI):  
                #     ZDat[i]= savgol_filter(ZDat[i], 65, 5)
                for ii in range(3):       
                    ar0x=np.mean(ZDat,axis=0)
                    ar0x_=0.4*(np.mean(abs((ZDat)-(ar0x)),axis=0))
                        
                    lnn=len(ZDat[0])
                    NNew=int(.35*lnn)
                    for i in range(anI):    
                        for j in range(lnn):    
                            if not abs(ZDat[i,j]-(ar0x[j]))<=ar0x_[j]:     
                                if ZDat[i,j]<((ar0x[j])-ar0x_[j]):            
                                    ZDat[i,j]=(ar0x[j])-ar0x_[j]
                                else:
                                    if ZDat[i,j]>((ar0x[j])+ar0x_[j]):
                                        ZDat[i,j]=(ar0x[j])+ar0x_[j]

                    P=np.zeros(3,float)
                    for i in range(anI):
                        dd=ZDat[i][lnn-NNew:].copy()                         
                        x=ar0x[lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))].copy()
                        ZDat[i][lnn-NNew:]=filterFourierQ(ZDat[i],(ar0x),NNew,1)[lnn-NNew:]
                        # ZDat[i]= savgol_filter(ZDat[i], 65, 5)
                        P[0:2]=np.polyfit(x,ZDat[i][lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))],1)
                        if not P[0]>0:
                            P[0:2]=np.polyfit(dd,ZDat[i][lnn-NNew:],1)
                        ZDat[i][lnn-NNew:]=(ZDat[i][lnn-NNew:]-P[1])/P[0]    
                        # ZDat[i]= savgol_filter(ZDat[i], 14, 5)    
                        
                ar0x[Nf-NNew:]=savgol_filter(ar0x[Nf-NNew:], 14, 5)     
                
                if Lo:
                    ar0x=np.exp(ar0x) 
                ar0x[0:len(ar0)]=ar0[0:len(ar0)].copy()

                if Lo:
                    dd_=np.log(ar0x)                   
                else:
                    dd_=ar0x.copy()
                    
                astart0=np.inf
                
                if Lo:
                    ar0_=np.exp(dd_)
                else:
                    ar0_=dd_.copy()
                MMM_=0
                tm0=tm.time()
                tm1=0
                aMx0=dd_*0-np.inf
                aMn0=dd_*0+np.inf
                mm1=dd_*0
                mm2=mm1.copy()
                while MMM_<2*Nproc and (tm1-tm0)<MxTime:              
                    arr_RezM=  np.zeros((Ngroup,Nf),float)  
                    arr_RezN=  np.zeros((Ngroup,Nf),float)  
                    MMM=0    
                    for iGr in range(Ngroup):   
                        ZDat=(Arr_AAA[iGr*NIter*int(Nproc/Ngroup)+((hh0))*int(Nproc/Ngroup):iGr*NIter*int(Nproc/Ngroup)+(hh0+min(dNIt, hhha-hh0))*int(Nproc/Ngroup)]).copy()
                        if Lo:
                            ZDat=np.exp(ZDat)
                        xxxx=ZDat.copy()
                        for i in range(aNN-1):
                            xxxx=np.concatenate((xxxx, ZDat))
                        ZDat=xxxx.copy()
                        hhhx=0
                        anI=len(ZDat)
                        # for i in range(anI):  
                        # #     ZDat[i][:len(ar0)]=ar0.copy()
                        #     if Lo:
                        #         ZDat[i]= np.exp(savgol_filter(np.log(ZDat[i]), 34, 5))
                        #     else:
                        #         ZDat[i]= savgol_filter(ZDat[i], 34, 5)
                        lnn=len(ZDat[0])
                        for i in range(anI):    
                            for j in range(lnn):    
                                if Lo:     
                                    ZDat[i]=np.log(ZDat[i])
                                    if not abs(ZDat[i,j]-np.log(ar0x[j]))<=ar0x_[j]:     
                                        if ZDat[i,j]<(np.log(ar0x[j])-ar0x_[j]):            
                                            ZDat[i,j]=np.log(ar0x[j])-ar0x_[j]
                                        else: 
                                            if ZDat[i,j]>(np.log(ar0x[j])+ar0x_[j]):
                                                ZDat[i,j]=np.log(ar0x[j])+ar0x_[j]
                                    ZDat[i]=np.exp(ZDat[i])
                                else:
                                    if not abs(ZDat[i,j]-(ar0x[j]))<=ar0x_[j]:     
                                        if ZDat[i,j]<((ar0x[j])-ar0x_[j]):            
                                            ZDat[i,j]=(ar0x[j])-ar0x_[j]
                                        else:
                                            if ZDat[i,j]>((ar0x[j])+ar0x_[j]):
                                                ZDat[i,j]=(ar0x[j])+ar0x_[j]
                                
                        if hhh==0:
                            if Lo:
                                arr_rezBzz=np.exp(np.mean(np.log(ZDat),axis=0)) 
                            else:
                                arr_rezBzz=np.mean(ZDat,axis=0)
                        for i in range(anI):
                            if not astart0==np.inf:
                                if Lo:                                
                                    ZDat[i][1:]=np.diff(np.log(ZDat[i])+KPP*np.log(arr_rezBzz))/(1+KPP)
                                else:
                                    ZDat[i][1:]=np.diff(ZDat[i]+KPP*arr_rezBzz)/(1+KPP)
                                ZDat[i][0]=0     
                            else:
                                if Lo:
                                    ZDat[i]=(np.log(ZDat[i])+KPP*np.log(arr_rezBzz))/(1+KPP)
                                else:
                                    ZDat[i]=(ZDat[i]+KPP*arr_rezBzz)/(1+KPP)
                        P=np.zeros(3,float)
                        for i in range(anI):
                            dd=ZDat[i][Nf-NNew:].copy()                         
                            if Lo:
                                x=np.log(ar0_[Nf-NNew:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))])
                                ZDat[i][Nf-NNew:]=filterFourierQ(ZDat[i],np.log(ar0_),NNew,1,-1,1)[Nf-NNew:]
                            else:
                                x=ar0_[Nf-NNew:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))].copy()
                                ZDat[i][Nf-NNew:]=filterFourierQ(ZDat[i],(ar0_),NNew,1,-1,1)[Nf-NNew:]
                            P[0:2]=np.polyfit(x,ZDat[i][Nf-NNew:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))],1)
                            if not P[0]>0:
                                # P[0:2]=np.polyfit(dd,ZDat[i][Nf-NNew:],1)   
                                P[0]=np.std(ZDat[i][Nf-NNew:])/np.std(dd)
                                P[1]=np.mean(ZDat[i][Nf-NNew:])-P[0]*np.mean(dd)
                            else:
                                P[0]=np.std(ZDat[i][Nf-NNew:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))])/np.std(x)
                                P[1]=np.mean(ZDat[i][Nf-NNew:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))])-P[0]*np.mean(x)   
                            
                            ZDat[i][Nf-NNew:]=savgol_filter((ZDat[i][Nf-NNew:]-P[1])/P[0], 14, 5)     
                        
                        if anI<aNN: 
                            all_RezM[iGr][hhh]=np.amax(ZDat,axis=0)
                            all_RezN[iGr][hhh]=np.amin(ZDat,axis=0)
                        else:
                            aMx_=0
                            aMn_=0
                            mdd4=ZDat*0
                            aa=RandomQ(Nf)                        
                            ss0=np.concatenate((aa, aa, aa))
                            hhhx=0
                            #pool = mp.Pool(processes=Numproc)
                            while hhhx<int(NIter/25):
                                #if dNIt*int(hhhx/dNIt)==hhhx:
                                areza=[]        
                                
                                ssssss=np.asarray([Nf,anI,NNew,hhhx,aMM,aNN,Numproc],np.int32)
                                aSred=np.mean(ZDat)
                                sssss=np.asarray(ZDat.reshape((anI*Nf))-aSred,np.float32)
                                ssss=np.asarray(ss0,np.int32)
                                
                                arrgsx=(struct.pack("%si"%(len(ssssss)),*ssssss)+
                                            struct.pack("%sf"%(anI*Nf),*sssss)+
                                            struct.pack("%si"%(len(ss0)),*ssss)) 
                                with open("maxmindat.tmp", "wb") as file:
                                    file.write(arrgsx)
                                    file.close()          
                                # RALF1Cella(zip(arrgsx))                                  
                                
                                wwwww=0
                                while not wwwww: 
                                    with open("maxmin.tmp", "wb") as file:
                                        file.close()  
                                    pool = mp.Pool(processes=Numproc)
                                    resa=pool.map_async(RALF1Cella,zip(arrgsx))
                                    wwwwww=0
                                    while not wwwwww:
                                        with open("maxmin.tmp", "rb") as file:                    
                                            areza=file.read()
                                            file.close()    
                                            areza=np.asarray(struct.unpack('%sf'%(int(len(areza)/4)), areza),np.float32)
                                        try:
                                            areza=np.asarray(areza.reshape((int(len(areza)/2/Nf),2,Nf)),np.float32)+aSred
                                            nsz=len(areza)
                                            
                                            if nsz>=Numproc:
                                                wwwww=1
                                                wwwwww=1 
                                                pool.terminate()
                                            else:
                                                tm.sleep(0.1)                                        
                                        except:
                                            wwwwww=1 
                                            pool.terminate()

                                #del(pool)

                                areza=np.asarray(areza,float)
                                PP=1
                                aMx=np.amax(areza[:,0],axis=0)
                                aMn=np.amin(areza[:,1],axis=0)
                                
                                aMx=savgol_filter(aMx, 14, 5)
                                aMn=savgol_filter(aMn, 14, 5)
                                
                                if Lo:
                                    x=np.log(ar0_[:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))])
                                else:
                                    x=ar0_[:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))].copy()
                                    
                                if Lo:
                                    y_1=aMx[:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))].copy()
                                    y_2=aMn[:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))].copy()

                                else:
                                    y_1=aMx[:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))].copy()
                                    y_2=aMn[:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))].copy()
                              
                                P_1=P.copy()
                                P_2=P.copy()
                                try:
                                    P_1[0:2]=np.polyfit(x,y_1,1)
                                    P_2[0:2]=np.polyfit(x,y_2,1)
                                    
                                    # P_1[0]=np.std(y_1)/np.std(x)
                                    # P_1[1]=np.mean(y_1)-P_1[0]*np.mean(x)
                                    # P_2[0]=np.std(y_2)/np.std(x)
                                    # P_2[1]=np.mean(y_2)-P_2[0]*np.mean(x)

                                    if not P_1[0]>0 or not P_2[0]>0 or 100*scp.pearsonr(x,y_1)[0]<10 or 100*scp.pearsonr(x,y_2)[0]<10:
                                        PP=0                                        
                                except:
                                    PP=0
                                
                                if not PP==0: 
                                    aMx=(aMx-P_1[1])/P_1[0]
                                    aMn=(aMn-P_2[1])/P_2[0]  
                                    if dNIt*int(hhhx/dNIt)==hhhx:
                                        aMx_=aMx.copy()
                                        aMn_=aMn.copy()
                                        aMx0=aMx_.copy()
                                        aMn0=aMn_.copy()
                                        
                                    aMx_=(aMx_*(hhhx-dNIt*int(hhhx/dNIt))+np.maximum(aMx_,aMx))/(hhhx-dNIt*int(hhhx/dNIt)+1)
                                    aMn_=(aMn_*(hhhx-dNIt*int(hhhx/dNIt))+np.minimum(aMn_,aMn))/(hhhx-dNIt*int(hhhx/dNIt)+1)
                                    aMx0=aMx_.copy()
                                    aMn0=aMn_.copy()
                                                                                               
                                    if not PP==0:
                                        if hhhx==0:
                                            dd1=aMx_.copy()
                                            dd2=aMn_.copy()
                                        dd1=np.maximum(dd1,aMx_)
                                        dd2=np.minimum(dd2,aMn_)
                                        hhhx=hhhx+1
                                    else:
                                        PP=0
                                if PP==0:
                                    aMx_=aMx0.copy()
                                    aMn_=aMn0.copy()
                                    aa=RandomQ(Nf)                        
                                    ss0=np.concatenate((aa, aa, aa))
                                # tm1=tm.time()
                                # if (tm1-tm0)>MxTime:
                                #     break
                            del(pool)  
                            
                        # tm1=tm.time()
                        # if (tm1-tm0)>MxTime:
                        #     break  
                                              
                        if Lo:
                            x=np.log(ar0_[:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))])
                        else:
                            x=ar0_[:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))].copy()
                        
                        y_1=dd1[:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))].copy()
                        y_2=dd2[:Nf-NNew+int(lSrez*(NNew-(Nf-len(ar0))))].copy()
                        
                        P_1[0:2]=np.polyfit(x,y_1,1)
                        P_2[0:2]=np.polyfit(x,y_2,1)
                        
                        # P_1[0]=np.std(y_1)/np.std(x)
                        # P_1[1]=np.mean(y_1)-P_1[0]*np.mean(x)
                        # P_2[0]=np.std(y_2)/np.std(x)
                        # P_2[1]=np.mean(y_2)-P_2[0]*np.mean(x)
                        
                        arr_RezM[iGr][Nf-NNew:]=0.5*((dd1[Nf-NNew:]-P_1[1])/P_1[0]+
                                                           (dd2[Nf-NNew:]-P_2[1])/P_2[0])
                      
                        if Lo:
                            arr_RezM[iGr][:Nf-NNew]=np.log(ar0[:Nf-NNew])
                        else:
                            arr_RezM[iGr][:Nf-NNew]=ar0[:Nf-NNew].copy()
                            
                        arr_RezM[iGr]= savgol_filter(arr_RezM[iGr], 14, 5)
                        all_RezMM[iGr][hhh]=arr_RezM[iGr].copy()                        
                        arr_RezM[iGr]=np.mean(all_RezMM[iGr][:hhh+1],axis=0)
                                                    
                    # tm1=tm.time()
                    # if (tm1-tm0)>MxTime:
                    #     break 
                    MMM=int(2*MMM/Ngroup)
                    arr_rezBz=np.mean(arr_RezM, axis=0) 
 
                    if Lo:
                        arr_rezBz=np.exp(arr_rezBz) 
                        for iGr in range(Ngroup): 
                            arr_RezM[iGr]=np.exp((arr_RezM[iGr])) 
                            
                    mm1=ar0[Nf-NNew:].copy()                            
                    mm2=arr_rezBz[Nf-NNew:len(ar0)].copy()   
                    # try: 
                    #     if 100*scp.pearsonr(mm1,mm2)[0]>10 and MMM==0:
                    #         break
                    #     else:
                    #         MMM_=MMM_+1
                    # except:
                    #     break
                    break
                
                # tm1=tm.time()
                # if MMM_<2*Nproc and np.std(mm1)>0 and np.std(mm2)>0 and (tm1-tm0)<MxTime:
                kcoef=100*scp.pearsonr(mm1,mm2)[0]
                
                try:
                    if (kcoef>0.7*Koef_[len(Koef_)-1]):
                        Koef_.append(kcoef) 
                    else:
                        MMM_=2*Nproc
                except:
                        Koef_.append(kcoef) 
                                    
                if MMM_<2*Nproc:# and np.std(mm1)>0 and np.std(mm2)>0 and (tm1-tm0)<MxTime:
                    arr_rezBzz=arr_rezBz.copy()
                    
                    
                    y_data = arr_rezBzz.copy()
                    y_data[:len(arrr)]=arrr.copy()
                    x_data = np.asarray(range(len(y_data)))
                    mask = (x_data >=len(arrr)-1) * (x_data<int(Nf*.98))
                    x_range = x_data[mask]
                    y_range = y_data[mask]
                    
                    max_idx = np.argmax(y_range)
                    try:
                        min_idx = np.argmin(y_range[:max_idx])
                        
                        min_y = y_range[min_idx]
                        max_y = y_range[max_idx]
                        min_x = x_range[min_idx]
                        max_x = x_range[max_idx]
                        if min_x < max_x:
                            profit_percentage1 = (max_y - min_y) / abs(min_y) * 100
                        else:
                            profit_percentage1 = np.nan
                    except:
                        profit_percentage1 = np.nan
                        
                    min_idx = np.argmin(y_range)
                    try:
                        max_idx = np.argmax(y_range[min_idx:])
                        
                        min_y = y_range[min_idx]
                        max_y = y_range[max_idx]
                        min_x = x_range[min_idx]
                        max_x = x_range[max_idx]
                        if min_x > max_x:
                            profit_percentage2 = (max_y - min_y) / abs(min_y) * 100
                        else:
                            profit_percentage2 = np.nan
                    except:
                        profit_percentage2 = np.nan
                        
                    profit_percentage=max(profit_percentage1, profit_percentage2)

                    fig = plt.figure()
                    axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
                    axes_ = fig.add_axes([0, 0, 0.3, 0.3]) 
                    axes__ = fig.add_axes([0.4, 0, 0.3, 0.3])        
                    
                    if Lo:
                        axes.semilogy(ar0, 'ro-', alpha=0.1)
                        axes.semilogy(arrr, 'rx-')
                        for iGr in range(Ngroup):
                            axes.semilogy(arr_RezM[iGr,:int(Nf*.98)],linewidth=3.,alpha=0.2)
                        axes.semilogy(arr_rezBz[:int(Nf*.98)],'yx-',linewidth=4.,alpha=0.5)
                        axes.grid(True, which="both", ls="-")
                    else:
                        axes.plot(ar0, 'ro-', alpha=0.1)
                        axes.plot(arrr, 'rx-')
                        for iGr in range(Ngroup):
                            axes.plot(arr_RezM[iGr,:int(Nf*.98)],linewidth=3.,alpha=0.2)
                        axes.plot(arr_rezBz[:int(Nf*.98)],'yx-',linewidth=4.,alpha=0.5)
                        axes.grid(True, which="both", ls="-")
                    
                    axes.text(-0.1, 4, '%s  '%(hhh+1),
                            verticalalignment='bottom', horizontalalignment='left',
                            transform=axes_.transAxes,color='black', fontsize=24) 
                    axes.text(4, 4, 'Reflection on Action of Lorentz Forces-1, #2011612714  \n\n Course = %s, start = %s, step = %s * %s / %s'%(aname,adat0,interv,aDecm,andecim),
                            verticalalignment='bottom', horizontalalignment='right',
                            transform=axes_.transAxes,color='blue', fontsize=14)    
                    if Lo:
                        gint=np.polyfit(np.log(mm1),np.log(mm2),1)
                    else:
                        gint=np.polyfit(mm1,mm2,1)
                    
                    if Lo:                
                        axes_.loglog(mm1,np.exp(gint[1]+gint[0]*np.log(mm1)),'y.',linewidth=2.)                    
                        axes_.loglog(mm1,mm2, 'ok', markersize=3, alpha=0.1) 
                    else:
                        axes_.plot(mm1,gint[1]+gint[0]*mm1,'y.',linewidth=2.)                    
                        axes_.plot(mm1,mm2, 'ok', markersize=3, alpha=0.1) 
                    
                    axes_.text(0.2, 0.6, '%d'%int(np.floor(np.asarray(Koef_,float)[::-1][0])),
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes_.transAxes,color='green', fontsize=30)
                    
                    axes__.text(1.8, 0.6, 'Dunning-Kruger\n effect',
                            verticalalignment='bottom', horizontalalignment='center',
                        transform=axes_.transAxes,color='green', fontsize=14)  
                    axes__.plot(np.asarray(range(hhh+1),float)+1,Koef_,'y',linewidth=2.)
                    
                    try:
                        # Add percentage in bottom right corner
                        axes.text(4.2, 0.4, '%d%%'%int(np.floor(0.7*profit_percentage)),
                            verticalalignment='bottom', horizontalalignment='right',
                            transform=axes_.transAxes,color='magenta', fontsize=50)
                    except:
                        # Add percentage in bottom right corner
                        axes.text(4, -0.2, '-',
                            verticalalignment='bottom', horizontalalignment='right',
                            transform=axes_.transAxes,color='magenta', fontsize=160)

                        
                    frame=fig2img(fig) 
                    cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)        
                    gray_sz1=min(gray_sz1,len(cimg[0]))
                    gray_sz2=min(gray_sz2,len(cimg))
                    ImApp.append(frame)
                    #or 3*int((hhh+1)/5)==(hhh+1) or 
                    if WrtTodr>0 or hhha==(hhh+1):
                        out = cv.VideoWriter(wrkdir + aname+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
                        for icl in range(len(ImApp)):
                            cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
                            out.write(cimgx[0:gray_sz2,0:gray_sz1,:]) 
                        out.release()
                    del(out)
                    plt.show()
                    hhh=hhh+1

                else:
                    try:
                        dill.load_session(wrkdir + aname+".ralf")
                    except:
                        hh0=hh0    
                hh0=hh0+1
                #if WrtTodr>0:
                try:
                    dill.dump_session(wrkdir + aname+".ralf")  
                except:
                    hh0=hh0
                if hh0==2*NIter:
                    hhh=NIter 
                print (hhh+10000*hh0)

# import pandas as pd
# import numpy as np
# import csv
# import scipy.interpolate as interp

# Lo = 1
# # Original WhO list with duplicate ICP (handled in processing)
# WhO = [
#     'SHBBTC','XRP', 'ryujin', "AAVE","ANKR","AVAX","AXS","BAL","BTC","BCH","ADA","CGLD","LINK","COMP",
#     "CRV","ICP","MANA","DOGE","EOS","ETH","ETC","GRT","LTC","LRC","MKR","NU",
#     "OMG","ORN","MATIC","SHIB","SKL","SOL","XLM","SNX","XTZ","UNI"
# ]

# WhO= [
#     "SNX",
#     "LTC",
#     "BAT",
#     "BCH",
#     "ETC",
#     "ALGO",
#     "CRV",
#     "MANA",
#     "ATOM",
#     "DOGE",
#     "UNI",
#     "1INCH",
#     "ANKR",
#     "XTZ",
#     "LINK",
#     "AVAX",
#     "BTC",
#     "GRT",
#     "ADA",
#     "XLM",
#     "MATIC",
#     "DOT",
#     "XRP",
#     "SHIB",
#     "AAVE",
#     "ICP",
#     "SOL",
#     "ETH",
#     "CELO",
#     "APE",
#     "CHZ"
# ]

# WhO= ['ANKR', 'BAT', 'MANA', 'CHZ', 'UNI', 'ETC', 'DOT', 'ATOM', 'ETH']
# wrkdir = r""
# GetSCV = 0
# aDecm = 1
# andecim=1
# Lengt0=300

# def decimat(adat_):
#     if Lo:
#         adat_ = np.log(adat_)
#     k = 0
#     adat__ = np.zeros(int(len(adat_) / aDecm), float)
#     for i in range(int(len(adat_) / aDecm)):
#         adat__[k] = np.median(adat_[i * aDecm:i * aDecm + aDecm])
#         k = k + 1

#     arrray=adat__[1:len(adat__)].copy()
    
#     indices = np.arange(len(arrray))
#     non_nan_mask = ~np.isnan(arrray)
#     nan_mask = np.isnan(arrray)
    
#     interpolated_values = np.interp(indices[nan_mask], indices[non_nan_mask], arrray[non_nan_mask])
    
#     arrray[nan_mask] = interpolated_values
    
#     if Lo:
#         return np.exp(arrray)
#     else:
#         return (arrray)
    
# def getcsv(WhO, xYears, wrkdir):
#     for i in range(len(WhO)):
#         nm = WhO[i]
#         anm = "https://www.cryptodatadownload.com/cdd/Binance_%sUSDT_d.csv" % (WhO[i])
#         try:
#             df = pd.read_csv(anm, delimiter=",", skiprows=[0])
#             arrr = df['Close']
#             arrr.to_csv(wrkdir + nm + '.csv', index=False)
#         except:
#             anm = anm
# if GetSCV:
#     getcsv(WhO, 5, wrkdir)

# def loaddata(aLengt, ticker1, key):
#     adat_ = []
#     with open(wrkdir + ticker1 + '.csv', newline='') as csvfile:
#         spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
#         dat = []
#         i = 0
#         for row in spamreader:
#             if i > 0:
#                 dat.append(row[0])
#             i = i + 1
#     dat = dat[len(dat)::-1]
#     try:
#         dat = np.asarray(dat, float)
#         arrr = np.asarray(dat, float)
#     except:
#         dat = dat
#         arrr = np.zeros(2, int)
        
#     if andecim>=1:
#         sz=len(arrr)
#         xold=np.asarray(range(sz),int)
#         f=interp.interp1d(xold, arrr,'cubic')
#         xnew=np.asarray(range((sz-1)*andecim),int)/andecim
#         arrr=f(xnew)
#         sz=len(arrr)
#         arrr=arrr[sz-aLengt*aDecm:sz]
#         # %varexp --plot arrr 
        
#     return arrr, adat_

# arrrxxR = []
# nams = []
# lenar = []

# for i in range(len(WhO)):
#     i0 = 0
#     while i0 < 6 and not i0 < 0:
#         try:
#             arrrxx1, adat1_ = loaddata(Lengt0, WhO[i], 1)
#             arrrxxR.append(arrrxx1)
#             lenar.append(len(arrrxx1))
#             nams.append(WhO[i])
#             i0 = -1
#         except:
#             i0 = i0 + 1
#             #tm.sleep(10)

# llar = int(np.median(np.asarray(lenar, int)))
# nnams_ = []
# aaer = []
# for i in range(len(nams)):
#     if len(arrrxxR[i]) >= llar:
#         aer = decimat(arrrxxR[i])
#         aaer.append(aer[len(aer) - int(llar / aDecm) + 1:].copy())
# aaer = np.asarray(aaer, float)

# arrrxxR_ = []
# ii = 0
# for i in range(len(nams)):
#     if len(arrrxxR[i]) >= llar:
#         # arrrxxR_.append(np.diff(np.log(aaer[ii])))
#         zzzz=np.log(aaer[ii])
#         zzzz=(zzzz-np.mean(zzzz))/np.std(zzzz)
#         zzzz_=zzzz
#         zzzz_[0]=0
#         zzzz_[1:]=np.diff(zzzz)
#         arrrxxR_.append(zzzz+zzzz_)        
#         nnams_.append(nams[ii])
#         ii = ii + 1
# arrrxxR_ = np.asarray(arrrxxR_, float)

# srarr = np.median((arrrxxR_), axis=0)
# arrrxxRR_ = np.asarray(arrrxxR_, float) * 0
# for i in range(len(nnams_)):
#     arrrxxRR_[i] = (arrrxxR_[i] - srarr)

# # Compute the correlation matrix
# correlation_mat = np.corrcoef(arrrxxRR_)

# # Sort based on average absolute correlation
# avg_correlations = np.mean(np.abs(correlation_mat), axis=0)
# sorted_indices = np.argsort(avg_correlations)  # Sort ascending (lowest to highest correlation)
# correlation_mat_sorted = correlation_mat[sorted_indices][:, sorted_indices]
# nnams_sorted = [nnams_[i] for i in sorted_indices]

# # Extract the most independent cryptocurrencies (top 25% lowest average correlation)
# num_independent = max(1, len(nnams_sorted) // 4)  # At least 1, up to 25% of the list
# most_independent_indices = sorted_indices[:num_independent]
# most_independent_coins = [nnams_sorted[i] for i in range(num_independent)]
# print("Most independent cryptocurrencies:", most_independent_coins)

# # Visualize the sorted correlation matrix
# import seaborn as sns
# import matplotlib.pyplot as plt
# sns.set(font_scale=2)
# plt.figure(figsize=(12, 10))
# sns.heatmap(correlation_mat_sorted, 
#             xticklabels=nnams_sorted,
#             yticklabels=nnams_sorted,
#             cmap='RdPu',
#             vmin=-0.2, vmax=1.0,
#             cbar_kws={'label': 'Correlation'})
# plt.title('Correlation Matrix of Cryptocurrencies (Sorted by Average Correlation)\n')
# plt.show()
################################
# import hickle as hkl
# import numpy as np
# from RALf1FiltrVID import filterFourierQ
# from scipy.signal import savgol_filter
# import pylab as plt
# from PIL import Image  
# import cv2 as cv
# import os
# def fig2dsk ( fig ):
#     fig.savefig(wrkdir +'_tm_dynamic.png',dpi=150,transparent=False,bbox_inches = 'tight')
#     frame=Image.open(wrkdir +'_tm_dynamic.png')
#     return frame
# lSrez=0.99
# arrrxx=hkl.load("ralfrez.rlf2")
# # i=0 
# # arrrxx_=[]
# # for ii in arrrxx:
# #     if not i==4:
# #         arrrxx_.append(ii)
# #     i=i+1
# # hkl.dump(arrrxx_,"ralfrez.rlf2")
# try:
#     ZDat=np.asarray(arrrxx,float)#.transpose()
#     anI=len(ZDat)
#     # for i in range(anI):  
#     #     ZDat[i]= savgol_filter(ZDat[i],65, 5)
#     for ii in range(3):       
#         ar0x=np.median(ZDat,axis=0)
#         ar0x_=.4*(np.median(abs((ZDat)-(ar0x)),axis=0))
            
#         lnn=len(ZDat[0])
#         NNew=int(.35*lnn)
#         for i in range(anI):    
#             for j in range(lnn):    
#                 if not abs(ZDat[i,j]-(ar0x[j]))<=ar0x_[j]:     
#                     if ZDat[i,j]<((ar0x[j])-ar0x_[j]):            
#                         ZDat[i,j]=(ar0x[j])-ar0x_[j]
#                     else:
#                         if ZDat[i,j]>((ar0x[j])+ar0x_[j]):
#                             ZDat[i,j]=(ar0x[j])+ar0x_[j]
        
#         P=np.zeros(3,float)
#         for i in range(anI):
#             dd=ZDat[i][lnn-NNew:].copy()                         
#             x=ar0x[lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))].copy()
#             ZDat[i][lnn-NNew:]=filterFourierQ(ZDat[i],(ar0x),NNew,1)[lnn-NNew:]
#             # ZDat[i]= savgol_filter(ZDat[i],65, 5)
#             P[0:2]=np.polyfit(x,ZDat[i][lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))],1)
#             if not P[0]>0:
#                 P[0:2]=np.polyfit(dd,ZDat[i][lnn-NNew:],1)
#             ZDat[i][lnn-NNew:]=(ZDat[i][lnn-NNew:]-P[1])/P[0]
#             ZDat[i]= savgol_filter(ZDat[i], 14, 5)                      
#     bbbbb=ZDat[:,:].transpose().copy()
#     aaaaa=np.median(bbbbb.transpose(),axis=0)
        
#     %varexp --plot bbbbb 
#     len(arrrxx)
# except:
#     len(arrrxx)
# len(arrrxx)
# aname="UNI"

# fourcc = 0
# ImApp=[]
# gray_sz1=0

# dNIt=4
# Ngroup=3
# Lo=1
# Nproc=Ngroup*(os.cpu_count())
# wrkdir = r""
# [hhhao,Arr_AAA]=(hkl.load(wrkdir + aname+".rlf1"))
# NIter=100

# for hhhai in range(hhhao):  
#     hhha=hhhai+1
#     for iGr in range(Ngroup):  
#         ZDat=Arr_AAA[iGr*NIter*int(Nproc/Ngroup)+max(0,(hhha)-dNIt)*int(Nproc/Ngroup):iGr*NIter*int(Nproc/Ngroup)+(hhha)*int(Nproc/Ngroup)].copy()
#         if iGr==0:
#             xxxx=ZDat.copy()
#         else:
#             xxxx=np.concatenate((xxxx, ZDat))
#     ZDat=xxxx.copy()   
#     anI=len(ZDat)
#     # for i in range(anI):  
#     #     ZDat[i]= savgol_filter(ZDat[i], 65, 5)
#     for ii in range(3):       
#         ar0x=np.median(ZDat,axis=0)
#         ar0x_=0.4*(np.median(abs((ZDat)-(ar0x)),axis=0))
            
#         lnn=len(ZDat[0])
#         NNew=int(.35*lnn)
#         for i in range(anI):    
#             for j in range(lnn):    
#                 if not abs(ZDat[i,j]-(ar0x[j]))<=ar0x_[j]:     
#                     if ZDat[i,j]<((ar0x[j])-ar0x_[j]):            
#                         ZDat[i,j]=(ar0x[j])-ar0x_[j]
#                     else:
#                         if ZDat[i,j]>((ar0x[j])+ar0x_[j]):
#                             ZDat[i,j]=(ar0x[j])+ar0x_[j]

#         P=np.zeros(3,float)
#         for i in range(anI):
#             dd=ZDat[i][lnn-NNew:].copy()                         
#             x=ar0x[lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))].copy()
#             ZDat[i][lnn-NNew:]=filterFourierQ(ZDat[i],(ar0x),NNew,1)[lnn-NNew:]
#             # ZDat[i]= savgol_filter(ZDat[i], 65, 5)
#             P[0:2]=np.polyfit(x,ZDat[i][lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))],1)
#             if not P[0]>0:
#                 P[0:2]=np.polyfit(dd,ZDat[i][lnn-NNew:],1)
#             ZDat[i][lnn-NNew:]=(ZDat[i][lnn-NNew:]-P[1])/P[0]    
#             ZDat[i]= savgol_filter(ZDat[i], 14, 5)                  
#     #ZDat=Arr_AAA[iGr*NIter*int(Nproc/Ngroup)+max(0,(hhh+1)-dNIt)*int(Nproc/Ngroup):iGr*NIter*int(Nproc/Ngroup)+(hhh+1)*int(Nproc/Ngroup)].copy()
#     bbbbb=ZDat[:,:].transpose().copy()
#     aaaaa=(np.amax(bbbbb.transpose(),axis=0)+np.amin(bbbbb.transpose(),axis=0))/2
#     fig = plt.figure()
#     axes = fig.add_axes([0.2, 0.2, 0.8, 0.8])    
#     axes.plot(bbbbb)
#     frame=fig2dsk(fig) 
#     cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)     
#     if gray_sz1==0:       
#         gray_sz1=len(cimg[0])
#         gray_sz2=len(cimg)
#     else:
#         gray_sz1=min(gray_sz1,len(cimg[0]))
#         gray_sz2=min(gray_sz2,len(cimg))
#     ImApp.append(frame)
#     if fourcc==0:
#         aDur=3
#         fourcc = cv.VideoWriter_fourcc(*'MP4V')
#     out = cv.VideoWriter(wrkdir + aname+'_tm_'+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
#     for icl in range(len(ImApp)):
#         cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
#         out.write(cimgx[0:gray_sz2,0:gray_sz1,:]) 
#     out.release()
#     del(out)
#     plt.show()

