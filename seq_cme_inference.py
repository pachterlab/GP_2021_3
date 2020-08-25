
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
# %matplotlib inline

import time
import loompy as lp
import scipy
import velocyto as vcy
# %config InlineBackend.figure_format = 'retina'

import time
import os
from datetime import date
import pickle
from multiprocessing import Pool
import random

import sklearn
from sklearn.cluster import KMeans
import collections
import warnings
	

########################
## Statistical testing
########################
def chisq_gen(result_data,viz=False,nosamp=False):
    NGEN_TESTED = result_data.n_gen
    import scipy.stats.mstats
    if not nosamp:
        expected_freq = [cme_integrator_samp(
                    np.insert(10**result_data.best_phys_params[i_],0,1), 
                    10**result_data.gene_spec_samp_params[i_], 
                    result_data.M[i_],result_data.N[i_],np.inf,'none').flatten() for i_ in range(NGEN_TESTED)]
    else:
        expected_freq = [cme_integrator_nosamp(
                    np.insert(10**result_data.nosamp_gene_params[i_],0,1), 
                    result_data.M[i_],result_data.N[i_],np.inf,'none').flatten() for i_ in range(NGEN_TESTED)]
    for i_ in range(NGEN_TESTED):
        temp = expected_freq[i_]
        EPS=1e-12
        temp[temp<EPS]=EPS
        expected_freq[i_] = temp
    csqarr = [scipy.stats.mstats.chisquare(result_data.hist[i_].flatten(), 
                                           expected_freq[i_]) for i_ in range(NGEN_TESTED)]
    csq = np.array([csqarr[i_][0] for i_ in range(len(csqarr))])
    pval = np.array([csqarr[i_][1] for i_ in range(len(csqarr))])

    result_data.set_rej(pval,nosamp=nosamp)

    if viz:
        plt.hist(csq)
        plt.xlabel('Chi-square statistic')
        plt.ylabel('# genes')
    return (csq,pval)

########################
## Analysis
########################

def import_datasets(pickle_filenames):
    #result .pickle file(s), must be given as list
    result_data = ResultData()
    for i_ in range(len(pickle_filenames)):
        with open(pickle_filenames[i_],'rb') as hf:
            SO = pickle.load(hf)
            if i_==0: #note this presupposes all are defined over same grid, potentially different genes
                result_data.set_parameters(SO)
            result_data.set_variables(SO)
    return result_data

def landscape_viz(RES,log=True,colorbar=True):
    sz = (RES.n_pt2,RES.n_pt1)
    if log:
        DIVG = np.log10(RES.divg)
    else:
        DIVG = RES.divg
    
    X = np.reshape(RES.X.T,sz)
    Y = np.reshape(RES.Y.T,sz)
    Z = np.reshape(DIVG.T,sz)
    plt.contourf(X,Y,Z,40)
    if colorbar:
        plt.colorbar()

def resample_opt_viz(RES,resamp_vec=(1,2,3,4,5),Ntries=4,figsize=(10,10)):
    sz = (RES.n_pt2,RES.n_pt1)
    Nsamp = len(resamp_vec)

    fig1,ax1=plt.subplots(nrows=Nsamp,ncols=Ntries,figsize=figsize)
    for samp_num in range(Nsamp):
        for i_ in range(Ntries):
            axis_location = (samp_num,i_) #np.unravel_index(plot_index,sz)

            gene_selection = np.random.choice(RES.n_gen,resamp_vec[samp_num],replace=False)

            divg_samp = np.sum(RES.gene_spec_err[:,gene_selection],1)
            
            X = np.reshape(RES.X.T,sz)
            Y = np.reshape(RES.Y.T,sz)
            Z = np.reshape(divg_samp.T,sz)
            
            ax1[axis_location].contourf(X,Y,np.log10(Z),15)
            loc_best_ind = np.argmin(divg_samp)
            ax1[axis_location].scatter(RES.X[loc_best_ind],RES.Y[loc_best_ind],s=10,c='r')
            
            
            if i_==0:
                ax1[axis_location].set_ylabel('Ngen = '+str(resamp_vec[samp_num]))

            ax1[axis_location].set_xticks([])
            ax1[axis_location].set_yticks([])

def resample_opt_mc_viz(RES,resamp_vec=(1,2,3,4),Ntries=1000,figsize=(8,2),log=True):
    sz = (RES.n_pt2,RES.n_pt1)
    
    if log:
        DIVG = np.log10(RES.divg)
    else:
        DIVG = RES.divg
    Nsamp = len(resamp_vec)
    
    fig1,ax1=plt.subplots(nrows=1,ncols=Nsamp,figsize=figsize)
    # for plot_index in range(N_):
    for samp_num in range(Nsamp):
        LOC = []
        for i__ in range(Ntries):

            gene_selection = np.random.choice(RES.n_gen,resamp_vec[samp_num],replace=False)

            divg_samp = np.sum(RES.gene_spec_err[:,gene_selection],1)
            loc_best_ind = np.argmin(divg_samp)
            LOC.append([RES.X[loc_best_ind],RES.Y[loc_best_ind]])       
        LOC = np.asarray(LOC)
    
        X = np.reshape(RES.X.T,sz)
        Y = np.reshape(RES.Y.T,sz)
        Z = np.reshape(DIVG.T,sz)

        ax1[samp_num].contourf(X,Y,Z,30)
        jit = np.random.normal(scale=0.03,size=LOC.shape)
        LOC=LOC+jit
        ax1[samp_num].scatter(LOC[:,0],LOC[:,1],c='r',s=3,alpha=0.3)
        ax1[samp_num].set_xticks([])
        ax1[samp_num].set_yticks([])
        ax1[samp_num].set_title('Ngen = '+str(resamp_vec[samp_num]))

def plot_param_marg(result_data,nbin=15,nosamp=False):
    from scipy.stats import norm

    fig1,ax1=plt.subplots(nrows=1,ncols=3,figsize=(5,2))

    param_nm = ('burst size','deg rate','splice rate')
    for i in range(3):
        if not nosamp:
            DATA = result_data.best_phys_params[:,i]
            LB = result_data.phys_lb[i]
            UB = result_data.phys_ub[i]
        else: 
            DATA = result_data.nosamp_gene_params[:,i]
            LB = result_data.phys_lb_nosamp[i]
            UB = result_data.phys_ub_nosamp[i]

        ax1[i].hist(DATA,nbin,density=True)
    #     print(np.mean(best_phys_params[:,i]))

        mu, std = norm.fit(DATA)
        
        xmin, xmax = ax1[i].get_xlim()
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, std)
        ax1[i].plot(x, p, 'k', linewidth=2)
        
        ax1[i].set_xlim([LB,UB])
        ax1[i].set_title(param_nm[i])
        ax1[i].set_xlabel('log10 value')
    fig1.tight_layout()

def plot_param_L_dep(result_data,nosamp=False):
    fig1,ax1=plt.subplots(nrows=1,ncols=3,figsize=(5,2))

    name_var = ('log b','log gamma','log beta')
    for i in range(3):
        if not nosamp:
            DATA = result_data.best_phys_params[:,i]
            LB = result_data.phys_lb[i]
            UB = result_data.phys_ub[i]
        else: 
            DATA = result_data.nosamp_gene_params[:,i]
            LB = result_data.phys_lb_nosamp[i]
            UB = result_data.phys_ub_nosamp[i]

        ax1[i].scatter(result_data.gene_log_lengths,DATA,c='k',s=1,alpha=0.5)
        ax1[i].set_xlabel('log L')
        ax1[i].set_ylabel(name_var[i])
        ax1[i].set_ylim([LB,UB])
    fig1.tight_layout()

def plot_KL(result_data,nbins=15,nosamp=False):
    if not nosamp:
    	DATA = result_data.gene_spec_err[result_data.best_ind]
    else:
    	DATA = result_data.nosamp_gene_spec_err
    plt.hist(DATA,nbins)

    plt.xlabel('KL divergence')
    plt.ylabel('# genes')

def plot_genes(result_data,sz,figsize,marg='none',log=False,title=True,nosamp=False):
    (nrows,ncols)=sz
    fig1,ax1=plt.subplots(nrows=nrows,ncols=ncols,figsize=figsize)
    
    NGEN_PLOT = np.prod(sz)
    for i_ in range(NGEN_PLOT):
        axis_location = np.unravel_index(i_,sz)
        
        if not nosamp:
            Pa = cme_integrator_samp(numpy.insert(10**result_data.best_phys_params[i_],0,1),
                                    10**result_data.gene_spec_samp_params[i_],
                                     result_data.M[i_],result_data.N[i_],np.inf,marg)
        else:
            Pa = cme_integrator_nosamp(numpy.insert(10**result_data.nosamp_gene_params[i_],0,1),
                                     result_data.M[i_],result_data.N[i_],np.inf,marg)

        if log and marg == 'none':
        	Pa[Pa<1e-8]=1e-8
        	Pa = np.log10(Pa)
        if title:
            if hasattr(result_data,'gene_rej') and result_data.gene_rej[i_] and not nonsamp:
                ax1[axis_location].set_title(result_data.gene_names[i_]+' (rejected)',fontdict={'fontsize': 9})
            if hasattr(result_data,'gene_rej_nosamp') and result_data.gene_rej_nosamp[i_] and nonsamp:
                ax1[axis_location].set_title(result_data.gene_names[i_]+' (rejected)',fontdict={'fontsize': 9})
            else:
                ax1[axis_location].set_title(result_data.gene_names[i_],fontdict={'fontsize': 9})
        ax1[axis_location].set_xticks([])
        ax1[axis_location].set_yticks([])
        if marg=='none':
            
            X_,Y_ = np.meshgrid(np.arange(result_data.M[i_])-0.5,
                                np.arange(result_data.N[i_])-0.5)
            ax1[axis_location].contourf(X_.T,Y_.T,Pa,20,cmap='summer')
            
            jitter_magn = 0.1
            jitter_x = np.random.randn(result_data.Ncells)*jitter_magn
            jitter_y = np.random.randn(result_data.Ncells)*jitter_magn
            ax1[axis_location].scatter(result_data.raw_U[i_]+jitter_x,
                                       result_data.raw_S[i_]+jitter_y,c='k',s=1,alpha=0.1)
            
            ax1[axis_location].set_xlim([-0.5,result_data.M[i_]-1.5])
            ax1[axis_location].set_ylim([-0.5,result_data.N[i_]-1.5])
        if marg=='nascent':
            ax1[axis_location].hist(result_data.raw_U[i_],
                                    bins=np.arange(result_data.M[i_])-0.5,density=True,log=log)
            ax1[axis_location].plot(np.arange(result_data.M[i_]),Pa)
            ax1[axis_location].set_xlim([-0.5,result_data.M[i_]-1.5])
        if marg=='mature':
            ax1[axis_location].hist(result_data.raw_S[i_],
                                    bins=np.arange(result_data.N[i_])-0.5,density=True,log=log)
            ax1[axis_location].plot(np.arange(result_data.N[i_]),Pa.T)
            ax1[axis_location].set_xlim([-0.5,result_data.N[i_]-1.5])


    fig1.tight_layout(pad=0.02)

def chisq_best_param_correction(result_data,method='nearest',Niter_=10,viz=True,szfig=(2,5),figsize=(10,3),overwrite=False):
    if viz:
        fig1,ax1=plt.subplots(nrows=szfig[0],ncols=szfig[1],figsize=figsize)

    divg_orig = result_data.divg
    best_params = np.zeros((Niter_,2))
    for i_ in range(Niter_):
        (chisq,pval) = chisq_gen(result_data)
        result_data.divg = np.sum(result_data.gene_spec_err[:,~result_data.gene_rej],1)
        result_data.find_best_params()

        best_params[i_,:] = result_data.best_samp_params

        if viz:
            axl = np.unravel_index(i_,szfig)
            sz = (result_data.n_pt2,result_data.n_pt1)
            X = np.reshape(result_data.X.T,sz)
            Y = np.reshape(result_data.Y.T,sz)
            Z = np.log10(np.reshape(result_data.divg.T,sz))
            ax1[axl].contourf(X,Y,Z,40)
            ax1[axl].scatter(result_data.best_samp_params[0],result_data.best_samp_params[1],s=10,c='r')

            ax1[axl].set_xticks([])
            ax1[axl].set_yticks([])

    #return everything to original values
    if not overwrite:
	    result_data.divg = divg_orig
	    result_data.find_best_params()
	    (chisq,pval) = chisq_gen(result_data)

    best_param_est = np.mean(best_params,0)
    if method == 'nearest':
        return result_data.sampl_vals[np.argmin(np.sum((np.array(result_data.sampl_vals)-best_param_est)**2,1))]
    if method == 'raw':
        return best_param_est


########################
## Initialization
########################

def create_dir(search_data,dataset_dir,ID,DATESTRING=date.today().strftime("%y%m%d"),code_ver='',model='BSD'):
    file_string = dataset_dir+'gg_'+DATESTRING+'_'+model+'_'+str(search_data.n_pt1)+'x'+str(search_data.n_pt2)+'_'+str(search_data.n_gen)+'gen_'+str(ID)
    search_data.set_file_string(file_string)
    try: 
        os.mkdir(file_string) 
        with open(file_string+'/metadata.pickle','wb') as hf:
        	pickle.dump(search_data,hf)
    
        print('Directory ' + file_string+ ' created; metadata written.')
    except OSError as error: 
        print('Directory ' + file_string+ ' exists.')
    return file_string

def build_grid(n_pts,samp_lb,samp_ub):
    X,Y = np.meshgrid(np.linspace(samp_lb[0],samp_ub[0],n_pts[0]),np.linspace(samp_lb[1],samp_ub[1],n_pts[1]))
    X=X.flatten()
    Y=Y.flatten()
    sampl_vals = list(zip(X,Y))
    return (X,Y,sampl_vals)

def get_transcriptome(transcriptome_filepath,repeat_thr=15):
    repeat_dict = {}
    len_dict = {}

    #remember that the convention is name - length - 5mers - 6mers -.... 50mers - more
    thr_ind = repeat_thr-3
    with open(transcriptome_filepath,'r') as file:   
        for line in file.readlines():
            d = [i for i in line.split(' ') if i]
            repeat_dict[d[0]] =  int(d[thr_ind])
            len_dict[d[0]] =  int(d[1])
    return (len_dict,repeat_dict)
    
def select_gene_set(loom_filepaths,feat_dict,viz=False,
                          results_to_exclude=[],seed=6,n_gen=10,
                          filt_param=(0.01,0.01,350,350,4,4)):

    n_datasets = len(loom_filepaths)

    for i_data in range(n_datasets):
        loom_filepath = loom_filepaths[i_data]
        print('Dataset: '+loom_filepath)
        vlm = vcy.VelocytoLoom(loom_filepath)
        gene_names_vlm = vlm.ra['Gene']
        #check which genes I have length data for
        sel_ind_annot = [k for k in range(len(gene_names_vlm)) if gene_names_vlm[k] in feat_dict]
        
        NAMES = [gene_names_vlm[k] for k in range(len(sel_ind_annot))]
        COUNTS = collections.Counter(NAMES)
        sel_ind = [x for x in sel_ind_annot if COUNTS[gene_names_vlm[x]]==1]

        print(str(len(gene_names_vlm))+' features observed, '+str(len(sel_ind_annot))+' match genome annotations. '
            +str(len(sel_ind))+' are unique. '
            +str(len(vlm.ca[list(vlm.ca.keys())[0]]))+' cells detected.')

        Ncells = len(vlm.ca[list(vlm.ca.keys())[0]])


        gene_names = gene_names_vlm[sel_ind]
        S_max = np.amax(vlm.S[sel_ind,:],1)
        U_max = np.amax(vlm.U[sel_ind,:],1)
        S_mean = np.mean(vlm.S[sel_ind,:],1)
        U_mean = np.mean(vlm.U[sel_ind,:],1)
        
        #aesthetics
        sz=(12,4)
        alf=0.15
        ptsz=3

        len_arr = np.asarray([feat_dict[k] for k in gene_names])

        warnings.filterwarnings("ignore")
        clusters = KMeans(n_clusters=2,init=np.asarray([[4,-2.5],[4.5,-0.5]]),algorithm="full").fit(
            np.vstack((np.log10(len_arr),np.log10(S_mean + 0.001))).T)
        warnings.resetwarnings()

        gene_cluster_labels = clusters.labels_
        
        if viz:
            var_name = ('S','U')
            var_arr = (S_mean,U_mean)

            fig1, ax1 = plt.subplots(nrows=1,ncols=2,figsize=sz)
            for i in range(2):
                ax1[i].scatter(np.log10(len_arr), np.log10(var_arr[i] + 0.001),s=ptsz,
                            c=gene_cluster_labels,alpha=alf,cmap="Spectral")
                ax1[i].set_xlabel('log10 gene length')
                ax1[i].set_ylabel('log10 (mean '+var_name[i]+' + 0.001)')
        
        gene_filter = np.array(gene_cluster_labels,dtype=bool)
        
        if viz:
            fig2, ax2 = plt.subplots(nrows=1,ncols=2,figsize=sz)
            for i in range(2):
                ax2[i].scatter(np.log10(len_arr)[gene_filter], 
                            np.log10(var_arr[i] + 0.001)[gene_filter],s=ptsz,c='k',alpha=alf)
                ax2[i].set_xlabel('log10 gene length')
                ax2[i].set_ylabel('log10 (mean '+var_name[i]+' + 0.001)')
                
        print(str(sum(gene_filter))+' genes retained as high-expression.')
        gene_filter2  = gene_filter \
            & (U_mean > filt_param[0]) \
            & (S_mean > filt_param[1]) \
            & (S_max < filt_param[2]) \
            & (U_max < filt_param[3]) \
            & (S_max > filt_param[4]) \
            & (U_max > filt_param[5])
        
                    
        gene_names_filt = gene_names[gene_filter2]
        vlm_gene_filter =  np.asarray([True if x in gene_names_filt else False for x in vlm.ra['Gene']],dtype=bool)
        vlm.filter_genes(by_custom_array=vlm_gene_filter)
        print(str(len(vlm.ra['Gene']))+' genes retained in loom structure based on filter.')
        
        random.seed(a=seed)

        sample_domain = np.arange(len(vlm.ra['Gene']))
        
        if len(results_to_exclude)>0:
            GN=[]
            for i_ in range(len(results_to_exclude)):
                with open(results_to_exclude[i_],'rb') as hf:
                    SO = pickle.load(hf)
                    GN.extend(SO.gene_names)
            print(str(len(GN))+' genes previously run...')
            GN = set(GN)
            print(str(len(GN))+' genes were unique.')
            sample_domain = [i_ for i_ in sample_domain if vlm.ra['Gene'][i_] not in GN]
        
        print(str(len(sample_domain))+' genes retained in loom structure based on previous results.')
        SAMPLE_DOMAIN_NAMES = vlm.ra['Gene'][sample_domain]
        if i_data == 0:
            set_intersection = set(SAMPLE_DOMAIN_NAMES)
        else:
            set_intersection = set_intersection.intersection(SAMPLE_DOMAIN_NAMES)
        print('Gene set size: '+str(len(set_intersection)))
        print('-----------')
        

    trunc_gene_set = np.array(list(set_intersection))
    if n_gen < len(trunc_gene_set):
        gene_select = np.random.choice(trunc_gene_set,n_gen,replace=False)
        print(str(n_gen)+' genes selected.')
    else:
        gene_select = trunc_gene_set
        print(str(len(trunc_gene_set))+' genes selected: cannot satisfy query of '+str(n_gen)+' genes.')
    
    gene_select=list(gene_select)
    trunc_gene_set = list(trunc_gene_set)
    return gene_select, trunc_gene_set

def get_gene_data(loom_filepath,feat_dict,gene_set,trunc_gene_set,viz=False):
    n_gen = len(gene_set)

    vlm = vcy.VelocytoLoom(loom_filepath)
    gene_names_vlm = vlm.ra['Gene']
    n_gen_tot = len(gene_names_vlm)
    #check which genes I have length data for
    sel_ind_annot = [k for k in range(len(gene_names_vlm)) if gene_names_vlm[k] in feat_dict]
    
    NAMES = [gene_names_vlm[k] for k in range(len(sel_ind_annot))]
    COUNTS = collections.Counter(NAMES)
    sel_ind = [x for x in sel_ind_annot if COUNTS[gene_names_vlm[x]]==1]
    Ncells = len(vlm.ca[list(vlm.ca.keys())[0]])

    print(str(len(gene_names_vlm))+' features observed, '+str(len(sel_ind_annot))+' match genome annotations. '
        +str(len(sel_ind))+' are unique. '
        +str(Ncells)+' cells detected.')


    sel_ind_bool_filter = np.zeros(n_gen_tot,dtype=bool)
    sel_ind_bool_filter[sel_ind] = True
    vlm.filter_genes(by_custom_array=sel_ind_bool_filter)
    gene_names = list(vlm.ra['Gene'])
    S_mean = np.mean(vlm.S,1)
    U_mean = np.mean(vlm.U,1)
    
    #aesthetics
    sz=(12,4)
    alf=0.15
    ptsz=3

    len_arr = np.asarray([feat_dict[k] for k in gene_names])

    warnings.filterwarnings("ignore")
    clusters = KMeans(n_clusters=2,init=np.asarray([[4,-2.5],[4.5,-0.5]]),algorithm="full").fit(
        np.vstack((np.log10(len_arr),np.log10(S_mean + 0.001))).T)
    warnings.resetwarnings()
    gene_cluster_labels = clusters.labels_
    
    gene_set_ind = [gene_names.index(gene_set[i_]) for i_ in range(n_gen)]
    trunc_gene_set_ind = [gene_names.index(trunc_gene_set[i_]) for i_ in range(len(trunc_gene_set))]
    low_expr_ind = np.where(gene_cluster_labels==0)[0]
    high_expr_filt_out = np.setdiff1d(
        np.where(gene_cluster_labels==1)[0],
        trunc_gene_set_ind)

    I_ = [low_expr_ind,high_expr_filt_out,trunc_gene_set_ind,gene_set_ind] 
    COL = [[0.9]*3, [0.8]*3, [0.2]*3, [0,0,1]]
    COL2 = [[0.3]*3, [0]*3, [0]*3]
    ALF = [0.2,0.3,0.3,0.8]
    ALF2 =  [0.02,0.1,0.1]


    if viz:
        
        warnings.filterwarnings("ignore")
        var_name = ('S','U')
        var_arr = (S_mean,U_mean)


        fig1, ax1 = plt.subplots(nrows=1,ncols=2,figsize=sz)
        for i in range(2):
            for j in range(3):
                ax1[i].scatter(np.log10(len_arr[I_[j]]), np.log10(var_arr[i][I_[j]] + 0.001),s=ptsz,
                            color=COL2[j],alpha=ALF2[j],cmap="Spectral")
            ax1[i].set_xlabel('log10 gene length')
            ax1[i].set_ylabel('log10 (mean '+var_name[i]+' + 0.001)')

        fig2, ax2 = plt.subplots(nrows=1,ncols=2,figsize=sz)
        for i in range(2):
            for j in range(4):
                ax2[i].scatter(np.log10(len_arr[I_[j]]), np.log10(var_arr[i][I_[j]] + 0.001),s=ptsz,
                            color=COL[j],alpha=ALF[j],cmap="Spectral")
            ax2[i].set_xlabel('log10 gene length')
            ax2[i].set_ylabel('log10 (mean '+var_name[i]+' + 0.001)')
        warnings.resetwarnings()
    
    gene_select = gene_set_ind
    M = np.asarray([int(np.amax(vlm.U[gene_index])+2) for gene_index in gene_select])
    N = np.asarray([int(np.amax(vlm.S[gene_index])+2) for gene_index in gene_select])

    hist = []
    moment_data = []
    gene_names = []
    gene_log_lengths = []
    raw_U = []
    raw_S = []
    for i_ in range(n_gen):
        H, xedges, yedges = np.histogram2d(vlm.U[gene_select[i_]],vlm.S[gene_select[i_]], 
                                          bins=[np.arange(M[i_]+1)-0.5,
                                          np.arange(N[i_]+1)-0.5],
                                          density=True)
        hist.append(H)

        #u var, u mean, s mean
        moments = [np.var(vlm.U[gene_select[i_]]), np.mean(vlm.U[gene_select[i_]]), np.mean(vlm.S[gene_select[i_]])]
        moment_data.append(moments)

        raw_U.append(vlm.U[gene_select[i_]])
        raw_S.append(vlm.S[gene_select[i_]])

        gene_names.append(vlm.ra['Gene'][gene_select[i_]])
        gene_log_lengths.append(np.log10(feat_dict[vlm.ra['Gene'][gene_select[i_]]]))
    raw_U = np.array(raw_U)
    raw_S = np.array(raw_S)
    gene_log_lengths = np.array(gene_log_lengths)
    moment_data = np.array(moment_data)
    
    search_data = SearchData()
    search_data.set_gene_data(M,N,hist,moment_data,gene_log_lengths,n_gen,gene_names,Ncells,raw_U,raw_S)
    
    return search_data
    
def dump_results(file_string,include_nosamp=False):
    divg = []
    T_ = []
    gene_params = []
    gene_spec_err = []

    with open(file_string+'/metadata.pickle','rb') as hf:
        SO = pickle.load(hf)
    for i in range(SO.N_pts):
        with open(file_string+'/grid_point_'+str(i)+'.pickle','rb') as hf:
            PK = pickle.load(hf)
            divg.append(PK[0])
            T_.append(PK[1])
            gene_params.append(PK[2])
            gene_spec_err.append(PK[3])
    divg = np.array(divg)
    T_ = np.array(T_)
    gene_params = np.array(gene_params)
    gene_spec_err = np.array(gene_spec_err)
    SO.set_results(divg,T_,gene_params,gene_spec_err)

    if include_nosamp:
	    with open(file_string+'/nosamp.pickle','rb') as hf:
	    	PL = pickle.load(hf)
	    	SO.set_nosamp_results(PL[2],PL[3])

    with open(file_string+'/result.pickle','wb') as hf:
    	pickle.dump(SO, hf)

########################
## Estimation for a given sequencing parameter tuple
########################
def nonvec_driver(search_data,i):
    SAMP_ = search_data.sampl_vals[i]
    ZZ = kl_obj(search_data,SAMP_)
    

    with open(search_data.file_string+'/grid_point_'+str(i)+'.pickle','wb') as hf:
	    pickle.dump((ZZ[0],ZZ[1],ZZ[2],ZZ[3],SAMP_,ZZ[4],
	        ('Obj func total','Runtime','Best transcriptional parameters','Obj func separate','Sample value','Init and final time')), 
	        hf)

def kl_obj(search_data,log_samp_fit_params):
    time_in = time.time()

    gene_itr = range(search_data.n_gen)

    log_samp_params = np.asarray(
        [(search_data.gene_log_lengths[i_] + log_samp_fit_params[0], 
          log_samp_fit_params[1]) for i_ in range(search_data.n_gen)])
    interior_params = (search_data.interior_search_restarts,
                       search_data.phys_lb,search_data.phys_ub,
                       search_data.interior_maxiter,search_data.init_pattern)

    param_list = [(log_samp_params[i_], search_data.hist[i_], 
                   search_data.M[i_], search_data.N[i_], 
                   interior_params, search_data.moment_data[i_]) for i_ in gene_itr]

    results = [gene_specific_optimizer(param_list[i_],i_) for i_ in gene_itr]
    errors = [results[i_][0] for i_ in gene_itr]
    x_arr = [results[i_][1] for i_ in gene_itr]
    obj_func =  sum(errors)
    # print(np.str(np.round(log_samp_fit_params,2))+'\t'+str(np.round(obj_func,2)))
    
    time_out = time.time()
    d_time = time_out-time_in
    
    return (obj_func,d_time,x_arr,errors,(time_in,time_out))

def gene_specific_optimizer(params,gene_index):
#     print('gene optimizer')
    samp_params = 10**params[0]
    gene_specific_data = params[1]
    M = params[2]
    N = params[3]
    interior_params = params[4]
    GSMD = params[5] #gene specific moment data
    time_in = time.time()
    num_restarts = interior_params[0]
    lb_log = interior_params[1]
    ub_log = interior_params[2]
    maxiter = interior_params[3]
    init_pattern = interior_params[4]
    
    lb = 10**lb_log
    ub = 10**ub_log
    
    bnd = scipy.optimize.Bounds(lb_log,ub_log)
    
    #u var, u mean, s mean
    b_fit = np.clip((GSMD[0] / GSMD[1] - 1) / samp_params[0] - 1, lb[0], ub[0])
    gamma_fit = np.clip(b_fit * samp_params[1] / GSMD[2], lb[1], ub[1])
    beta_fit = np.clip(b_fit * samp_params[0] / GSMD[1], lb[2], ub[2])
    x0 = np.log10(np.asarray([b_fit, gamma_fit, beta_fit]))
    if init_pattern == 'random':
    	x0 = np.random.rand(3)*(ub_log-lb_log)+lb_log
    
    res_arr = scipy.optimize.minimize(lambda x: kl_div(
        gene_specific_data,
        cme_integrator_samp(np.insert(10**x,0,1), samp_params, M,N,np.inf,'none')),
                                            x0=x0, bounds=bnd,
                                      options={'maxiter':maxiter,'disp':False})
    x = res_arr.x
    err = res_arr.fun
    err_orig = err
    
    ERR_THRESH = 0.99
    
    for rest_ind in range(num_restarts-1):
        x0_rand = np.random.rand(3)*(ub_log-lb_log)+lb_log
#         print(str(rest_ind) +'\t'+str(x0_rand))
        res_arr = scipy.optimize.minimize(lambda x: kl_div(
            gene_specific_data,
            cme_integrator_samp(np.insert(10**x,0,1), samp_params, M,N,np.inf,'none')),
                                                x0=x0_rand, bounds=bnd,
                                          options={'maxiter':maxiter,'disp':False})
        if res_arr.fun < err*ERR_THRESH:
            x = res_arr.x
            err = res_arr.fun
    return (err,x)

def cme_integrator_samp(phys,samp,M,N,t,marg):
    if marg=='mature':
        M=1
    elif marg=='nascent':
        N=1
    NN = N
    N = 1+int(np.ceil((NN-1)/2))
    l=np.arange(M)
    k=np.arange(N)
    u_ = np.exp(samp[0]*(np.exp(-2j*np.pi*l/M)-1))-1
    v_ = np.exp(samp[1]*(np.exp(-2j*np.pi*k/NN)-1))-1
    u,v=np.meshgrid(u_,v_)
    u=u.flatten()
    v=v.flatten()

    
    fun = lambda x: INTFUNC_(x,phys[1:],u,v)

    T_ = time.time()
    INT = scipy.integrate.quad_vec(fun,0,t,epsabs=1e-20,epsrel=1e-5,full_output=True)
    T__ = time.time()-T_
    I2 = INT[0]
    
    I = np.exp(I2*phys[0])
    I = np.reshape(I.T,(N,M))
    
    N_fin = N+1 if np.mod(NN-1,2)==0 else -1
    I=np.vstack((I,np.hstack((np.reshape(
        np.flipud(I[1:N_fin,0].conj()),(N_fin-2,1)),np.flip(I[1:N_fin,1:].conj(),(0,1))))))
    return np.real(np.fft.ifft2(I)).T 

def kl_div(data, proposal):
    EPS = 1e-12; 
    proposal[proposal<EPS]=EPS
    filt = data>0
    data = data[filt]
    proposal = proposal[filt]
    d=data*np.log(data/proposal)
    return np.sum(d)

def INTFUNC_(x,params,U,V):
    f = params[2]/(params[2]-params[1])
    Ufun = f*V*np.exp(-params[1]*x) + (U-V*f)*np.exp(-params[2]*x)
    filt = ~np.isfinite(f)
    if np.any(filt): #params[1] = params[2]
    	Ufun[filt] = np.exp(-params[1]*x)*(U[filt] + params[1]*V[filt]*x)
    Ufun = params[0]*Ufun
    return Ufun/(1-Ufun)

########################
## Non-sampling optimizer
########################

def nosamp_driver(search_data):
    ZZ = kl_obj_nosamp(search_data)

    with open(search_data.file_string+'/nosamp'+'.pickle','wb') as hf:
	    pickle.dump((ZZ[0],ZZ[1],ZZ[2],ZZ[3],(np.nan,np.nan),ZZ[4],
	        ('Obj func total','Runtime','Best transcriptional parameters','Obj func separate','Sample value','Init and final time')), 
	        hf)

def kl_obj_nosamp(search_data):
    time_in = time.time()

    gene_itr = range(search_data.n_gen)

    interior_params = (search_data.interior_search_restarts_nosamp,
                       search_data.phys_lb_nosamp,search_data.phys_ub_nosamp,
                       search_data.interior_maxiter,search_data.init_pattern)
    param_list = [([], search_data.hist[i_], 
                   search_data.M[i_], search_data.N[i_], 
                   interior_params, search_data.moment_data[i_]) for i_ in gene_itr]

    results = [gene_specific_optimizer_nosamp(param_list[i_],i_) for i_ in gene_itr]
    errors = [results[i_][0] for i_ in gene_itr]
    x_arr = [results[i_][1] for i_ in gene_itr]
    obj_func =  sum(errors)
    # print(np.str(np.round(log_samp_fit_params,2))+'\t'+str(np.round(obj_func,2)))
    
    time_out = time.time()
    d_time = time_out-time_in
    
    return (obj_func,d_time,x_arr,errors,(time_in,time_out))

def gene_specific_optimizer_nosamp(params,gene_index):
#     print('gene optimizer, no sampling')
    gene_specific_data = params[1]
    M = params[2]
    N = params[3]
    interior_params = params[4]
    GSMD = params[5] #gene specific moment data
    time_in = time.time()
    num_restarts = interior_params[0]
    lb_log = interior_params[1]
    ub_log = interior_params[2]
    maxiter = interior_params[3]
    init_pattern = interior_params[4]
    
    lb = 10**lb_log
    ub = 10**ub_log
    
    bnd = scipy.optimize.Bounds(lb_log,ub_log)
    
    #u var, u mean, s mean
    b_fit = np.clip(GSMD[0] / GSMD[1] - 1, lb[0], ub[0])
    gamma_fit = np.clip(b_fit / GSMD[2], lb[1], ub[1])
    beta_fit = np.clip(b_fit / GSMD[1], lb[2], ub[2])
    x0 = np.log10(np.asarray([b_fit, gamma_fit, beta_fit]))
    if init_pattern == 'random':
    	x0 = np.random.rand(3)*(ub_log-lb_log)+lb_log
    
    res_arr = scipy.optimize.minimize(lambda x: kl_div(
        gene_specific_data,
        cme_integrator_nosamp(np.insert(10**x,0,1), M,N,np.inf,'none')),
                                            x0=x0, bounds=bnd,
                                      options={'maxiter':maxiter,'disp':False})
    x = res_arr.x
    err = res_arr.fun
    err_orig = err
    
    ERR_THRESH = 0.99
    
    for rest_ind in range(num_restarts-1):
        x0_rand = np.random.rand(3)*(ub_log-lb_log)+lb_log
        res_arr = scipy.optimize.minimize(lambda x: kl_div(
            gene_specific_data,
            cme_integrator_nosamp(np.insert(10**x,0,1), M,N,np.inf,'none')),
                                                x0=x0_rand, bounds=bnd,
                                          options={'maxiter':maxiter,'disp':False})
        if res_arr.fun < err*ERR_THRESH:
            x = res_arr.x
            err = res_arr.fun
    return (err,x)

def cme_integrator_nosamp(phys,M,N,t,marg):
    if marg=='mature':
        M=1
    elif marg=='nascent':
        N=1
    NN = N
    N = 1+int(np.ceil((NN-1)/2))
    l=np.arange(M)
    k=np.arange(N)
    u_ = np.exp(-2j*np.pi*l/M)-1
    v_ = np.exp(-2j*np.pi*k/NN)-1
    u,v=np.meshgrid(u_,v_)
    u=u.flatten()
    v=v.flatten()
    
    fun = lambda x: INTFUNC_(x,phys[1:],u,v)

    T_ = time.time()
    INT = scipy.integrate.quad_vec(fun,0,t,epsabs=1e-20,epsrel=1e-5,full_output=True)
    T__ = time.time()-T_
    I2 = INT[0]
    
    I = np.exp(I2*phys[0])
    I = np.reshape(I.T,(N,M))
    
    N_fin = N+1 if np.mod(NN-1,2)==0 else -1
    I=np.vstack((I,np.hstack((np.reshape(
        np.flipud(I[1:N_fin,0].conj()),(N_fin-2,1)),np.flip(I[1:N_fin,1:].conj(),(0,1))))))
    return np.real(np.fft.ifft2(I)).T 


########################
## Class definitions
########################
class SearchData:
    def __init__(self):
        pass
    def set_gene_data(self,M,N,hist,moment_data,gene_log_lengths,n_gen,gene_names,Ncells,raw_U,raw_S):
        self.M = M
        self.N = N
        self.hist = hist
        self.moment_data = moment_data
        self.gene_log_lengths = gene_log_lengths
        self.n_gen = n_gen
        self.gene_names = gene_names
        self.Ncells = Ncells
        self.raw_U = raw_U
        self.raw_S = raw_S
    def set_interior_search_params(self,interior_search_restarts,phys_lb,phys_ub,interior_maxiter,init_pattern ='moments',phys_lb_nosamp=None,phys_ub_nosamp=None,interior_search_restarts_nosamp=None):
        self.interior_search_restarts = interior_search_restarts
        self.phys_lb = phys_lb
        self.phys_ub = phys_ub
        self.interior_maxiter = interior_maxiter
        self.init_pattern = init_pattern
        if phys_lb_nosamp is None:
        	phys_lb_nosamp = phys_lb
        if phys_ub_nosamp is None:
        	phys_ub_nosamp = phys_ub
        if interior_search_restarts_nosamp is None:
        	interior_search_restarts_nosamp = interior_search_restarts
        self.phys_lb_nosamp = phys_lb_nosamp
        self.phys_ub_nosamp = phys_ub_nosamp
        self.interior_search_restarts_nosamp = interior_search_restarts_nosamp

    def set_scan_grid(self,n_pt1,n_pt2,samp_lb,samp_ub):
        self.n_pt1 = n_pt1
        self.n_pt2 = n_pt2
        self.N_pts = n_pt1*n_pt2
        (X,Y,sampl_vals) = build_grid((n_pt1,n_pt2),samp_lb,samp_ub)
        self.X = X
        self.Y = Y
        self.sampl_vals = sampl_vals
    def set_file_string(self,file_string):
        self.file_string = file_string
    def get_pts(self):
        point_list = [i for i in range(self.N_pts) if not os.path.isfile(self.file_string+'/grid_point_'+str(i)+'.pickle')]
        print(str(len(point_list)) + ' of '+str(self.N_pts)+' points to be evaluated.')
        self.point_list = point_list
    def set_results(self,divg,T_,gene_params,gene_spec_err):
        self.divg = divg
        self.T_ = T_
        self.gene_params = gene_params
        self.gene_spec_err = gene_spec_err
    def set_nosamp_results(self,nosamp_gene_params,nosamp_gene_spec_err):
        self.nosamp_gene_params = nosamp_gene_params
        self.nosamp_gene_spec_err = nosamp_gene_spec_err

class ResultData:
    def __init__(self):
        self.gene_names = []
        self.hist = []
        self.M = np.zeros(0,dtype=int)
        self.N = np.zeros(0,dtype=int)
        self.gene_log_lengths = np.zeros(0)
        self.moment_data = np.zeros((0,3))
        self.n_gen = 0
        
    def set_parameters(self,search_results):
        attrs = ('n_pt1','n_pt2','N_pts','X','Y','sampl_vals',
                 'interior_search_restarts','interior_maxiter',
                 'phys_ub','phys_lb','Ncells')
        for attr in attrs:
            setattr(self,attr,getattr(search_results,attr))
        if has attr(search_results,'init_pattern'):
            setattr(self,'init_pattern',getattr(search_results,'init_pattern'))
        else:
            setattr(self,'init_pattern','moments')
        if hasattr(search_results,'phys_ub_nosamp'):
            setattr(self,'phys_ub_nosamp',getattr(search_results,'phys_ub_nosamp'))
            setattr(self,'phys_lb_nosamp',getattr(search_results,'phys_lb_nosamp'))
        self.divg = np.zeros(self.N_pts)
        self.gene_params = np.zeros((self.N_pts,0,3))
        self.gene_spec_err = np.zeros((self.N_pts,0))
        self.raw_U = np.zeros((0,self.Ncells),dtype=int)
        self.raw_S = np.zeros((0,self.Ncells),dtype=int)
        self.T_ = np.zeros((self.N_pts,0))

        if hasattr(search_results,'nosamp_gene_params'):
        	self.nosamp_gene_params = np.zeros((0,3))
        	self.nosamp_gene_spec_err = np.zeros(0)

    def set_variables(self,search_results):
        self.divg += search_results.divg
        self.gene_params = np.concatenate((self.gene_params,search_results.gene_params),axis=1)
        self.gene_spec_err = np.concatenate((self.gene_spec_err,search_results.gene_spec_err),axis=1)
        self.N = np.concatenate((self.N,search_results.N),axis=0)
        self.M = np.concatenate((self.M,search_results.M),axis=0)
        self.gene_log_lengths = np.concatenate((self.gene_log_lengths,search_results.gene_log_lengths),axis=0)
        self.moment_data = np.concatenate((self.moment_data,search_results.moment_data),axis=0)
        self.raw_U = np.concatenate((self.raw_U,search_results.raw_U),axis=0)
        self.raw_S = np.concatenate((self.raw_S,search_results.raw_S),axis=0)
        self.T_ = np.concatenate((self.T_,np.reshape(search_results.T_,(self.N_pts,1))),axis=1)
        self.gene_names.extend(search_results.gene_names)
        self.hist.extend(search_results.hist)
        self.n_gen += search_results.n_gen

        if hasattr(search_results,'nosamp_gene_params'):
        	self.nosamp_gene_params = np.concatenate((self.nosamp_gene_params,search_results.nosamp_gene_params),axis=0)
        	self.nosamp_gene_spec_err = np.concatenate((self.nosamp_gene_spec_err,search_results.nosamp_gene_spec_err),axis=0)

    def find_best_params(self):
        self.best_ind = np.argmin(self.divg)
        self.best_samp_params = self.sampl_vals[self.best_ind]
        self.best_phys_params = self.gene_params[self.best_ind]
        self.gene_spec_samp_params = np.array([(self.gene_log_lengths[i_] + self.best_samp_params[0], 
          self.best_samp_params[1]) for i_ in range(self.n_gen)])
    
    def set_rej(self,pval,threshold=0.05,bonferroni=True,nosamp=False):
        if bonferroni:
            threshold=threshold/self.n_gen
        if not nosamp:
            self.gene_rej = pval<threshold
        else:
        	self.gene_rej_nosamp = pval<threshold

