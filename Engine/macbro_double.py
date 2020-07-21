
def macbro_double(w,s,hwhmlist,hwhmlist2,p,mu):
    
    nw = len(w) # points in spectrum
    dw_new = (w[-1] - w[0]) / (nw-1) #wavelength change per pixel
    
    hwhm_max = max(hwhmlist)
    dw_min = dw_new
    
    #Pad spectrum ends to minimize impact of Fourier ringing.
    nhalf = abs(int(1.5*3.3972872 *hwhm_max/dw_min))
    ng = 2 * nhalf + 1 # points in gaussian (odd!)
    npad = nhalf + 2 # pad pixels on each end
    wg = dw_min * (np.arange(ng,dtype=float) - (ng-1)/2.0) #wavelength scale of gaussian
    xg = (0.83255461 / hwhm_max) * wg  #convenient absisca
    xg2 = (0.83255461 / hwhm_max) * (wg-p)  #convenient absisca
    
    spad = np.concatenate((s[0]*np.ones(npad),s,np.ones(npad)*s[-1]))
    hwhmpad = np.concatenate((hwhmlist[0]*np.ones(npad),hwhmlist,np.ones(npad)*hwhmlist[-1]))
    hwhm2pad = np.concatenate((hwhmlist2[0]*np.ones(npad),hwhmlist2,np.ones(npad)*hwhmlist2[-1]))
    
    sout = np.zeros_like(spad)
    
    for t in range(ng/2,len(spad)-ng/2):
        hwhm_new = hwhmpad[t]
        hwhm_new2 = hwhm2pad[t]
        
        gpro = (dw_new/hwhm_new)*(hwhm_max/dw_min)*np.exp(-xg*xg*((hwhm_max/hwhm_new)*(dw_new/dw_min))**2)
        gpro = gpro / np.sum(gpro)

        gpro2 = (dw_new/hwhm_new2)*(hwhm_max/dw_min)*np.exp(-xg2*xg2*((hwhm_max/hwhm_new2)*(dw_new/dw_min))**2)
        gpro2 = mu*gpro2 / np.sum(gpro2)

        
        gprotot = gpro + gpro2
        gprotot /= np.sum(gprotot)

        sout[t] = np.sum(spad[t-ng/2:t+ng/2+1]*gprotot)

    sout = sout[npad:npad+nw] #trim to original data/length
    sout = sout * np.sum(s) / np.sum(sout)
    return sout #return broadened spectrum
