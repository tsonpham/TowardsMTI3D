import tensorflow as tf
import numpy as np
import scipy
from pathlib import Path
from obspy import Stream, Trace
from obspy.geodetics import gps2dist_azimuth
import multiprocessing
import pickle
from functools import partial

_COMPs = {'r': 'r', 'theta': 't', 'phi': 'p'}
_EMTs = {'M_theta_theta':'Mtt', 'M_phi_phi':'Mpp', 'M_r_r':'Mrr', 
         'M_theta_phi':'Mtp', 'M_theta_r':'Mtr', 'M_phi_r':'Mpr'}

## reference Greens' tensor, which is known
def load_ses3d_gtensor(dname,vred=None,t0=None,wlen=None,dt=1):
    datapath = Path(dname)
    st = Stream()
    # objstats = []
    for emt in _EMTs:
        for fname in sorted((datapath/'DATA'/emt).glob('*.[xyz]')):
            netw, stat = fname.name.split('_')[0].split('.')
            loc = _EMTs[emt]
            ## parsing file name for network and station codes
            with open(str(fname)) as fp: content = fp.readlines()
            ## parse header info
            comp = _COMPs[content[0].split()[0]]
            npts = int(content[1].split()[1])
            delta = float(content[2].split()[1])
            xr, yr, zr = np.array(content[4].split()[1::2]).astype(float)
            xs, ys, zs = np.array(content[6].split()[1::2]).astype(float)
            dist, az, baz = gps2dist_azimuth(90-xs, ys, 90-xr, yr)
            dist /= 1e3
            ## cread data trace
            tr = Trace(header={'delta':delta}, data=np.loadtxt(content[7:]).astype(float))
            ## begin and end times
            if vred is not None:
                tbeg = dist/vred - t0
                tend = tbeg + wlen
                tr.trim(tr.stats.starttime+tbeg,tr.stats.starttime+tend,nearest_sample=False)
                # if (tr.stats.endtime-tr.stats.starttime < wlen):
                #     raise ValueError('Too long window length required.')
            ## append to stream
            st.append(tr)
    ## resample data and convert data into array
    st.resample(1/dt)
    st.filter('bandpass',freqmin=0.02,freqmax=0.05,zerophase=False,corners=3)
    ## convert to array
    nt = st[0].stats.npts
    ne,nc = len(_EMTs),3
    ns = int(len(st)/(ne*nc))
    gtensor = np.array([tr.data for tr in st]).reshape((ne,ns,nc,nt))
    gtensor = np.moveaxis(gtensor, 0, 2)
    return gtensor

def load_ses3d_gtensor_bulk(dname,vred=None,t0=None,wlen=None,dt=1):
    gf_path = dname/'GF3D.pkl'
    if not gf_path.is_file():
        # gtensors = [load_ses3d_gtensor(str(dname),vred,t0,wlen) for dname in dnames]
        with multiprocessing.Pool(4) as pool: 
            func = partial(load_ses3d_gtensor,vred=vred,t0=t0,wlen=wlen)
            gtensors = pool.map(func, sorted(dname.glob('MODEL???')))
        gtensors = np.array(gtensors)
        with open(str(gf_path), 'wb') as fp:
            pickle.dump(gtensors, fp)
    else:
        with open(str(gf_path), 'rb') as fp:
            gtensors = pickle.load(fp)
    return gtensors

def calc_Cov_t(m6, gtensors, gtensor_ref=None):
    nm, _, _, _, _ = gtensors.shape
    ## ensemble of synthetic waveforms
    syn = m6 @ gtensors
    if gtensor_ref is None:
        dev = syn - np.mean(syn, axis=0)
    else:
        dev = syn - m6 @ gtensor_ref
    ## theoretical covariance matrices
    Cov_t = np.mean(np.einsum('msct,mscu->msctu', dev, dev), axis=0) / (nm - 1)
    return Cov_t

def calc_inv_det(Cov, return_iLow=False):
    ns, nc, _, _ = Cov.shape
    ## Cholesky decomposition of low triangular matrix
    lower = np.linalg.cholesky(Cov)
    ## calculation of inversion of determinant
    iCov = np.zeros_like(Cov)
    log_Cov_det = np.zeros(Cov.shape[:-2])
    if return_iLow: iLow = np.zeros_like(Cov)
    for s in range(ns):
        for c in range(nc):
            lower_inv = scipy.linalg.inv(lower[s, c], True)
            if return_iLow: iLow[s, c] = lower_inv
            iCov[s, c] = lower_inv.T @ lower_inv
            log_Cov_det[s, c] = 2*np.sum(np.log(np.diag(lower[s, c])))
    ##
    if return_iLow:
        return iCov, log_Cov_det, iLow
    else:
        return iCov, log_Cov_det

def LinSolve(gtensor_ref, iCov, obs): 
    TMP = gtensor_ref @ iCov
    LHS = np.mean(TMP @ np.transpose(gtensor_ref, axes=(0, 1, 3, 2)), axis=(0, 1))
    RHS = np.mean(np.einsum('...ij,...j->...i', TMP, obs), axis=(0, 1))
    lsol = np.linalg.inv(LHS) @ RHS
    return lsol

def IterSolve(gtensor_ref, obs, noise_std, m6, gtensors, use_Ct_prime=False, n_iters=5, return_chain=False):
    nt = gtensors.shape[-1]
    m6_chain = []
    ## Iterative solution
    m6_chain.append(m6)
    for ii in range(n_iters):
        # structural uncertainty covariance matrix
        Cov_t = calc_Cov_t(m6, gtensors, gtensor_ref if use_Ct_prime else None)
        # combined covariance matrix with uncorrelated data noise
        Cov = np.einsum('sc...,sc->sc...', Cov_t, noise_std**-2) + np.eye(nt)
        # covariance matrices' inversion and determinant
        iCov, _ = calc_inv_det(Cov)
        iCov = np.einsum('sc...,sc->sc...', iCov, noise_std**-2)
        # linear solver
        m6 = LinSolve(gtensor_ref, iCov, obs)
        m6_chain.append(m6)
    if return_chain:
        return np.array(m6_chain)
    else:
        return m6_chain[-1]

def GradSolve(gtensor_ref, obs, learning_rate=2e-1, epochs=200, 
            M_init=None, t_init=None, random_seed=None, M_sigma=1, T_sigma=1.5):
    ns, _, ne, nt = gtensor_ref.shape
    ## gather data to creat Gtensor
    Gtensor = np.array(gtensor_ref)
    Gmed = np.median(np.abs(Gtensor))
    Obs =     np.array(obs)
    Omed = np.median(np.abs(Obs))
    ## preserve factor to scale back MT solution after inversion
    scale_factor = Omed / Gmed

    ########### tensorflow block, which will be run on GPU if available
    ## Fourier transforms of Gtensor and Obs waveforms 
    Gtensor   = tf.constant(Gtensor/Gmed, tf.float64)
    Gtensor_f = tf.signal.rfft(Gtensor, tf.constant([2*nt]))
    Obs       = tf.constant(Obs/Omed, tf.float64)
    Obs_f     = tf.signal.rfft(Obs, tf.constant([2*nt]))
    ## declare an optimizer
    optimizer = tf.keras.optimizers.Adam(learning_rate=learning_rate)

    ## randomised initial solutions
    if random_seed is not None: np.random.seed(random_seed)
    if M_init is None: 
        M = tf.Variable(np.random.normal(0, M_sigma, (1, ne)))
    else:
        M = tf.Variable(M_init.reshape((1, ne)))
    if t_init is None:
        t = tf.Variable(np.random.normal(0, T_sigma, (ns, 1, 1)))
    else:
        t = tf.Variable(t_init.reshape((ns, 1, 1)))

    ## frequency vector
    omega = tf.ones(Obs_f.shape, tf.float64) * tf.constant(np.fft.rfftfreq(2*nt)*2*np.pi, tf.float64)
    
    ## space holder to record loss and variable evolution
    history_loss, history_M, history_t = [], [], []
    for _ in range(epochs):
        # Open a GradientTape to record the operations run
        # during the forward pass, which enables auto-differentiation.
        with tf.GradientTape() as tape:
            ## observed data shifted instaneously by a time t
            Obs_f_real = tf.cos(omega*t)*tf.math.real(Obs_f) + tf.sin(omega*t)*tf.math.imag(Obs_f)
            Obs_f_imag = tf.cos(omega*t)*tf.math.imag(Obs_f) - tf.sin(omega*t)*tf.math.real(Obs_f)
            ## different of prediction and shifted observation
            Diff_real = tf.squeeze(M @ tf.math.real(Gtensor_f)) - Obs_f_real
            Diff_imag = tf.squeeze(M @ tf.math.imag(Gtensor_f)) - Obs_f_imag
            ## mean squared error
            loss_value = tf.reduce_mean(tf.square(Diff_real) + tf.square(Diff_imag))
            
        # Use the gradient tape to automatically retrieve
        # the gradients of the trainable variables with respect to the loss.
        grads = tape.gradient(loss_value, [M, t])
        # Run one step of gradient descent by updating
        # the value of the variables to minimize the loss.
        optimizer.apply_gradients(zip(grads, [M, t]))
        # Save history for inspection
        history_loss.append(float(loss_value))
        history_M.append(np.squeeze(M.numpy())*scale_factor)
        history_t.append(np.squeeze(t.numpy()))

    ## return a dictionary recording history evoluiton of loss function and its variables
    idx = np.argmin(history_loss)
    return history_M[idx]# , history_t[idx]

def SolveSuite(gtensor_ref, obs, noise_std, gtensors, inversion_type=dict()):
    '''
    Running the suite of inversion methods for MT inversion.
    Parameters:
    -----------
    gtensor_ref : np.ndarray
        Reference Green's tensors (ns, nc, ne, nt)
    obs : np.ndarray
        Observed data (ns, nc, nt)
    noise_std : np.ndarray
        Noise standard deviation (ns, nc)
    gtensors : np.ndarray
        Ensemble of Green's tensors (nm, ns, nc, ne, nt)
    inversion_type : dict
        Dictionary of inversion methods to be run, which includes:
        - 'wu_strucerror': iterative solution with unknown structural uncertainty
        - 'w_strucerror': inverted solution with true structural covariance matrices
        - 'wo_strucerror': MT solution without structure uncertainty
        - 'w_timeshift': MT solution without structure uncertainty but allowing time shifts
        the dictionary values contain necessary parameters for each inversion method.
    Returns:
    --------
    results : dict
        Dictionary of inversion results.
    '''
    ns, nc, _, nt = gtensor_ref.shape
    results = dict()
    if 'wu_strucerror' in inversion_type:
        ## interative solution with unknown structural uncertainty
        kw = inversion_type['wu_strucerror']
        m6 = np.random.uniform(-1, 1, 6)
        tmp = IterSolve(gtensor_ref, obs, noise_std, m6, gtensors, False, **kw)
        results['wu_C1_strucerror'] = tmp
        tmp = IterSolve(gtensor_ref, obs, noise_std, m6, gtensors, True, **kw)
        results['wu_C2_strucerror'] = tmp
    if 'w_strucerror' in inversion_type:
        ## inverted solution with true structural covariance matrices
        m6_true = inversion_type['w_strucerror']['m6_true']
        tmp = IterSolve(gtensor_ref, obs, noise_std, m6_true, gtensors, False, n_iters=1)
        results['w_C1_strucerror'] = tmp
        tmp = IterSolve(gtensor_ref, obs, noise_std, m6_true, gtensors, True, n_iters=1)
        results['w_C2_strucerror'] = tmp
    if 'wo_strucerror' in inversion_type:
        ## MT solution without structure uncertainty
        iCov = np.zeros((ns, nc, nt, nt))
        for s in range(ns):
            for c in range(nc):
                iCov[s, c] = np.eye(nt) * noise_std[s, c]**-2
        results['wo_strucerror'] = LinSolve(gtensor_ref, iCov, obs)
    if 'w_timeshift' in inversion_type:
        ## MT solution without structure uncertainty but allowing time shifts
        kw = inversion_type['w_timeshift']
        if 'wu_strucerror' in results: 
            M_init = results['wu_strucerror'][-1]
        else:
            M_init = np.zeros(6)
        results['w_timeshift'] = GradSolve(gtensor_ref, obs, M_init=M_init,
                                           t_init=np.zeros(ns), **kw)
    return results

import json
from hashlib import md5
import time

test1_dname = Path('InversionAssumptionTests')
def InversionAssumptionTests(root_path, m6_true, noise_scale=0.01):
    '''
    Test for multiple inversion assumptions (see Sections 2.4 and 2.5)
    Parameters:
    -----------
    '''
    ## Inversion assumptions
    inversion_type = dict(wu_strucerror={'n_iters':10},
                          w_strucerror={'m6_true':m6_true},
                          wo_strucerror={},
                          w_timeshift={'learning_rate':2e-1, 'epochs':500})
    
    ## Save results to file - to include inversion_type into the hashkey of the filename
    test1_dname.mkdir(parents=True, exist_ok=True)
    fname = md5(json.dumps([str(root_path), list(m6_true)]).encode('utf-8'))
    fname = test1_dname/fname.hexdigest()
    if fname.exists(): return

    ## Reference models of the known velocity models
    gtensor_ref = load_ses3d_gtensor('ses3d_r07_b', vred=vred, t0=t0, wlen=wlen)
    gtensor_ref = gtensor_ref[mask, :, :, :-1]

    ## Green's functions of randomised Earth's structures
    gtensors = load_ses3d_gtensor_bulk(root_path, vred, t0, wlen)
    gtensors = gtensors[:,mask,:,:,:-1]
    nm = gtensors.shape[0]
    
    ## Generate synthetic seismograms with uncorrelated data noise
    obs = m6_true @ gtensors
    noise_std = np.std(obs, axis=-1) * noise_scale
    obs += np.einsum('msct,msc->msct', np.random.normal(0, 1, obs.shape), noise_std)

    ## Repeat the inversion for true Earth model index
    results = [SolveSuite(gtensor_ref, obs[_], noise_std[_], gtensors, inversion_type) for _ in range(nm)]
    with open(fname, 'wb') as f: pickle.dump(results, f)

test2_dname = Path('IntegrationTests')
def IntegrationTests(root_path, m6_true_array, noise_scale=0.01):
    '''
    Test for multiple inversion assumptions (see Sections 2.4 and 2.5)
    Parameters:
    -----------
    '''    
    ## Save results to file - to include inversion_type into the hashkey of the filename
    test2_dname.mkdir(parents=True, exist_ok=True)
    fname = md5(json.dumps(str(root_path)).encode('utf-8'))
    fname = test2_dname/fname.hexdigest()
    if fname.exists(): return

    ## Reference models of the known velocity models
    gtensor_ref = load_ses3d_gtensor('ses3d_r07_b', vred=vred, t0=t0, wlen=wlen)
    gtensor_ref = gtensor_ref[mask, :, :, :-1]

    ## Green's functions of randomised Earth's structures
    gtensors = load_ses3d_gtensor_bulk(root_path, vred, t0, wlen)
    gtensors = gtensors[:,mask,:,:,:-1]
    nm = gtensors.shape[0]

    ## Repeat the process for multiple true MT solutions
    results = dict()
    for i, m6_true in enumerate(m6_true_array):
        print ('MT count', i)
        ## Inversion assumptions
        inversion_types = dict(wu_strucerror={'n_iters':10}, w_strucerror={'m6_true':m6_true})
        ## Generate synthetic seismograms with uncorrelated data noise
        obs = m6_true @ gtensors
        noise_std = np.std(obs, axis=-1) * noise_scale
        obs += np.einsum('msct,msc->msct', np.random.normal(0, 1, obs.shape), noise_std)
        ## Repeat the inversion for true Earth model index
        tmp = [SolveSuite(gtensor_ref, obs[_], noise_std[_], gtensors, inversion_types) for _ in range(nm)]
        results[tuple(m6_true)] = tmp

    ## Save pickled output to file for later uses
    with open(fname, 'wb') as f: pickle.dump(results, f)

if __name__ == '__main__':
    mask = np.array([1,1,0,1,1,1,0,0,1,0,0,0],dtype=bool)
    vred = 3.5 # km/s
    t0 = 35 # seconds after origin
    wlen = 115 # length of seismograms in seconds
    
    start = time.time()
    perturb_types = [
        'ses3d_r07_b/PERTURB_MODELS/UNCORRELATED_3D/1',
        'ses3d_r07_b/PERTURB_MODELS/UNCORRELATED_3D/3',
        'ses3d_r07_b/PERTURB_MODELS/UNCORRELATED_3D/5', 
        'ses3d_r07_b/PERTURB_MODELS/GAUSSIAN_3D/1',
        'ses3d_r07_b/PERTURB_MODELS/GAUSSIAN_3D/2',
        'ses3d_r07_b/PERTURB_MODELS/GAUSSIAN_3D/3',
        'ses3d_r07_b/PERTURB_MODELS/EXPONENTIAL_3D/1',
        'ses3d_r07_b/PERTURB_MODELS/EXPONENTIAL_3D/2',
        'ses3d_r07_b/PERTURB_MODELS/EXPONENTIAL_3D/3'
    ]
    # [InversionAssumptionTests(Path(_), m6_true) for _ in perturb_types]
    np.random.seed(1); m6_true_array = np.array([np.random.uniform(-1, 1, 6) for _ in range(100)])
    with multiprocessing.Pool(9) as p:
        p.starmap(IntegrationTests, [(Path(_), m6_true_array) for _ in perturb_types])
    print(f'Run time: {time.time()-start:.1f} seconds')