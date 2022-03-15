
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import copy
import pickle
from math import *
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import seaborn as sns
import scipy.stats as stats
import time
import multiprocessing
from multiprocessing import Pool, freeze_support

# g: maximum growth rate
# K: Monod constant
def OneStep(t, u, p):
    N, R, gmax, K, x = p
    Ni, Rj = u[:N], u[N:]
    xmat = np.hstack([x.T, 1-x.T])
    cmat = np.vstack([Rj for i in range(N)])
    gmat = xmat*(gmax*(cmat/(K+cmat)))
#     print(gmat)
    return np.hstack([np.sum(gmat, axis=1)*Ni, -Ni@gmat])


def Dilute(p, p_dilute, u_end):
    N, R, gmax, K, x = p
    c0, d, t, n, fluc_thres, elim_thres = p_dilute
    output = u_end
    u_add = np.zeros(N+R)
    u_add[-R:] = np.array(c0)
    stop, count = 1, 1
    utransf = copy.deepcopy(u_end)
    while stop>0 and count<n: 
        u = copy.deepcopy(utransf)
        u = u/d + u_add
        N_alive = u[:N]>0
        N_temp = sum(N_alive)
        gmax_temp = gmax[N_alive]
        K_temp = K[N_alive]
        u_temp = np.hstack([u[:N][N_alive], u[N:]])
        x_temp = x[0:,N_alive]
        p_temp = N_temp, R, gmax_temp, K_temp, x_temp
        sol = solve_ivp(OneStep, t_span=[0, t], y0=u_temp, method="LSODA", args=[p_temp], rtol=1e-6)
        utransf_temp = sol.y[:, -1] * (sol.y[:, -1]>elim_thres*d)
        u[:N][N_alive] = utransf_temp[:N_temp]
        u[-R:] = utransf_temp[-R:]
        utransf = u
        output = np.vstack([output, utransf])
        stop = sum(abs((output[-1, :-R]-output[-2, :-R])) > fluc_thres*output[-1, :-R])
        count+=1
#     print("c", count)
    return output

def invade(profile, p_profile, p_dilute, p_invade):
    # load variables
    N, R, gmax, K, x = p_profile
    c0, D, t, n, fluc_thres, elim_thres = p_dilute
    p_dilute_once = c0, D, t, 2, fluc_thres, elim_thres
    gmax_inv, K_inv, x_inv = p_invade
    p_profile_temp = N+1, R, np.vstack([gmax, gmax_inv]), np.vstack([K, K_inv]), np.hstack([x, x_inv])
    profile_temp = Dilute(p_profile_temp, p_dilute_once, np.hstack([profile[:-R], [D*elim_thres], profile[-R:]]))
    if(profile_temp[-1][-(R+1)]<=0):
        return profile, p_profile
    else:
        profile_temp = Dilute(p_profile_temp, p_dilute, np.hstack([profile[:-R], [5*D*elim_thres], profile[-R:]]))
        return profile_temp[-1], p_profile_temp

def model(index, N_pool):
    
    gA1, gA2, KA1, KA2 = 1, 1, 4, 5
    gB1, gB2, KB1, KB2 = 0.55, 0.55, 0.05, 0.03

    g_diff = np.random.normal(0, 0.05, size = 2*N_pool)
    K_exp = np.random.normal(0, 0.3, size = 4*N_pool)
    gA1new = [gA1+i for i in g_diff[0:N_pool]]
    gA2new = [gA2+i for i in g_diff[0:N_pool]]
    gB1new = [gB1+i for i in g_diff[N_pool:2*N_pool]]
    gB2new = [gB2+i for i in g_diff[N_pool:2*N_pool]]
    KA1new = [KA1*(10**i) for i in K_exp[0:N_pool]]
    KA2new = [KA2*(10**i) for i in K_exp[N_pool:2*N_pool]]
    KB1new = [KB1*(10**i) for i in K_exp[2*N_pool:3*N_pool]]
    KB2new = [KB2*(10**i) for i in K_exp[3*N_pool:4*N_pool]]
    
    xA, xB = 0.2, 0.7
    sigma = 0.1
    xAdist = stats.truncnorm(-0.2, 0.8, loc=0, scale=sigma)
    xAnew = np.sort(xA+xAdist.rvs(N_pool))[::-1]
    xBdist = stats.truncnorm(-0.7, 0.3, loc=0, scale=sigma)
    xBnew = np.sort(xB+xBdist.rvs(N_pool))
    
    plt.scatter(xAnew, 1-xAnew, color = 'b', alpha = 0.2)
    plt.scatter(xBnew, 1-xBnew, color = 'r', alpha = 0.2)

    N, R = 2*N_pool, 2
    gmax = np.array([[gA1new[i], gA2new[i]] for i in range(N_pool)]+[[gB1new[i], gB2new[i]] for i in range(N_pool)])
    K = np.array([[KA1new[i], KA2new[i]] for i in range(N_pool)]+[[KB1new[i], KB2new[i]] for i in range(N_pool)])
    x = np.array([np.hstack([xAnew, xBnew])])
    c0, D, t, n, fluc_thres, elim_thres = [10, 10], 10, 24.0, 5000, 1e-5, 1.0e-6
    
    start=time.time()
    
    p_profile = 1, R, gmax[:1, :], K[:1, :], x[:, :1]
    p_all = N, R, gmax, K, x
    p_dilute = c0, D, t, n, fluc_thres, elim_thres
    profile = Dilute(p_profile, p_dilute, np.array([1]+[0 for i in range(R)]))[-1]
    
    survivors_0 = [0]
    survivors = []
    while(survivors_0 != survivors):
        survivors_0 = survivors
        for i in range(N):
            _, _, _, _, x_profile = p_profile
            if(x[0, i] not in x_profile[0]):
                p_invade = gmax[i:i+1], K[i:i+1], x[:, i:i+1]
                profile_new, p_profile_new = invade(profile, p_profile, p_dilute, p_invade)
                N_temp, R_temp, gmax_temp, K_temp, x_temp = p_profile_new
                survivor = profile_new[:-R]>0
                p_profile = sum(survivor), R_temp, gmax_temp[survivor, :], K_temp[survivor, :], x_temp[:, survivor]
                profile = np.hstack([profile_new[:-R][survivor], profile_new[-R:]])

        end=time.time()

        N_temp, R_temp, gmax_temp, K_temp, x_temp = p_profile

        xlist = list(xAnew)+list(xBnew)
        survivors = [xlist.index(i) for i in x_temp[0]]

        all_survivors=survivors
        all_params=[xAnew, xBnew], gmax_temp, K_temp
        outputFileName = 'data/2_res_new_results_xandK_separated_'+str(index)+'.pkl'
        pickle.dump((all_survivors, all_params),open(outputFileName, "wb"))

    return (end-start)

runs = 500

if __name__ == '__main__':
    freeze_support()
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=7)
    x1 = [i for i in range(runs)]
    x2 = [50 for i in range(runs)]
    tasks = [(index, N_pool) for index, N_pool in zip(x1,x2)]
    pool.starmap(model,tasks)
