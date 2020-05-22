
# These classes are all for convenient passing of variables between the main code, the telluric spectrum fitter, and the model function.

import numpy as np

class fitobjs:

    def __init__(self,s, x, u, a0contwave,continuum,watm_in,satm_in,mflux_in,mwave_in):
        self.s = s
        self.x = x
        self.u = u
        self.a0contwave = a0contwave
        self.continuum  = continuum
        self.watm_in    = watm_in
        self.satm_in    = satm_in
        self.mflux_in   = mflux_in
        self.mwave_in   = mwave_in

class inparams:

    def __init__(self,inpath,outpath,initvsini,vsinivary,plotfigs,initguesses,bvcs,tagsA,tagsB,nights,mwave,mflux):
        self.inpath      = inpath
        self.outpath     = outpath
        self.initvsini   = initvsini
        self.vsinivary   = vsinivary
        self.plotfigs    = plotfigs
        self.initguesses = initguesses
        self.bvcs        = bvcs
        self.tagsA       = tagsA
        self.tagsB       = tagsB
        self.bvcs        = bvcs
        self.mwave0      = mwave
        self.mflux0      = mflux
        self.nights      = nights
        self.ips_tightmount = {2: np.array([2.04037987, 2.44724324, 1.94710219, 2.25831301, 2.04870909,
                                   2.54918255]), 3: np.array([2.5835808 , 2.51815444, 2.54622661, 2.33003733, 2.21602917,
                                   1.6183255 ]), 4: np.array([2.60463925, 2.50289512, 2.18010336, 2.08725543, 2.0447105 ,
                                   2.68329089]), 5: np.array([2.62590674, 2.34995633, 2.26398339, 2.10716032, 2.16699697,
                                   2.60612536]), 6: np.array([2.58312563, 2.4476741 , 2.28593311, 2.18518033, 2.50510778,
                                   2.63379758])}
        self.ips_loosemount = {2: np.array([5.23186917, 5.24900976, 3.96609885, 3.85262363, 3.11525618,
                                   3.6394157 ]), 3: np.array([5.68122268, 4.94600199, 4.42893859, 3.97849734, 3.1649909 ,
                                   2.35128276]), 4: np.array([5.49237893, 5.02570037, 3.91499657, 3.44131085, 3.04847074,
                                   3.3096322 ]), 5: np.array([5.63882553, 4.60165159, 4.01001552, 3.40237606, 3.06376238,
                                   2.98584427]), 6: np.array([5.25760553, 4.56458826, 3.93430453, 3.30376016, 3.11607041,
                                   3.0286326 ])}
        self.ips_tightmountFunc = {2: np.array([ 2.25378271e-09, -6.44641469e-06,  5.56041474e-03,  7.72955295e-01]),
                                    3: np.array([-1.03836176e-09,  2.28983215e-06, -1.75341703e-03,  2.97931406e+00]),
                                    4: np.array([ 2.20182527e-09, -5.41154039e-06,  3.32290632e-03,  2.00177689e+00]),
                                    5: np.array([ 1.00979196e-09, -1.94770701e-06,  3.99506340e-04,  2.68272193e+00]),
                                    6: np.array([ 1.52820030e-10,  4.24369854e-07, -1.35635083e-03,  3.04970584e+00])}

        self.ips_loosemountFunc = {2: np.array([ 4.05939305e-09, -1.11231891e-05,  7.06911774e-03,  3.98215719e+00]),
                                    3: np.array([-1.38630177e-09,  3.80492542e-06, -5.48547631e-03,  7.28241681e+00]),
                                    4: np.array([ 2.62862802e-09, -6.25435468e-06,  1.75570160e-03,  5.62385973e+00]),
                                    5: np.array([-7.39811208e-11,  1.81851999e-06, -5.51876063e-03,  7.46185799e+00]),
                                    6: np.array([ 8.12015898e-10, -1.16758952e-06, -2.22011098e-03,  6.23343180e+00])}

        # From GJ281 analysis with initrv = 20.02, not 20.15...
#        self.chunkweights_loosemount = {2: np.array([0.12484669, 0.01187605, 0.55409838, 0.15976596, 0.05187606, 0.09753685]),
#                                        3: np.array([0.15766037, 0.18344571, 0.16991371, 0.07209045, 0.0654897 , 0.35140006]),
#                                        4: np.array([0.09405195, 0.15956649, 0.1323817 , 0.28264784, 0.12322565, 0.20812638]),
#                                        5: np.array([0.05062509, 0.02441812, 0.04149529, 0.58755857, 0.07947067, 0.21643226]),
#                                        6: np.array([0.05591804, 0.22737467, 0.14304127, 0.13048763, 0.24373394, 0.19944446])}

#        self.chunkweights_tightmount = {2: np.array([0.20320015, 0.01013282, 0.64901208, 0.0412501 , 0.03356023, 0.06284462]),
#                                        3: np.array([0.2568777 , 0.20905393, 0.18243515, 0.06326623, 0.13547085, 0.15289615]),
#                                        4: np.array([0.16384428, 0.34125251, 0.16715084, 0.10815973, 0.15738878, 0.06220387]),
#                                        5: np.array([0.26709441, 0.22042356, 0.11256351, 0.00471484, 0.39341003, 0.00179364]),
#                                        6: np.array([0.09604141, 0.3705629 , 0.30773417, 0.01641322, 0.16791284, 0.04133546])}

        # From GJ281 analysis with Nsplit=8, initrv = 20.18: update: 04022020
        self.chunkweights_loosemount = {2: np.array([0.13992335, 0.00548394, 0.37635348, 0.21547168, 0.03589399, 0.22687356]),
                                        3: np.array([0.20653672, 0.16309499, 0.06580472, 0.04768457, 0.31613504, 0.20074395]),
                                        4: np.array([0.06278375, 0.18963667, 0.05826427, 0.25202664, 0.27911588, 0.1581728 ]),
                                        5: np.array([0.05326762, 0.02296659, 0.03711518, 0.33631157, 0.3878665 , 0.16247254]),
                                        6: np.array([0.03406971, 0.19150636, 0.15882698, 0.09657096, 0.28659474, 0.23243126])}

        self.chunkweights_tightmount = {2: np.array([0.21525525, 0.01333442, 0.61514991, 0.05957497, 0.01979944, 0.07688601]),
                                        3: np.array([0.21843997, 0.26091404, 0.13681933, 0.06954126, 0.13812044, 0.17616496]),
                                        4: np.array([0.1526929 , 0.29268455, 0.11404597, 0.19192788, 0.19955444, 0.04909426]),
                                        5: np.array([0.13061788, 0.07999278, 0.055955  , 0.3347575 , 0.22324851, 0.17542833]),
                                        6: np.array([0.05020937, 0.30509215, 0.33481064, 0.08170708, 0.16561447, 0.06256628])}
        self.methodvariance_tight = np.array([0.00282905, 0.00090454, 0.00124326, 0.00210822, 0.00140632])
        self.methodvariance_loose = np.array([0.00883368, 0.00840557, 0.00534123, 0.00385086, 0.00219453])







class inparamsA0:

    def __init__(self,inpath,outpath,plotfigs,tags,nights,humids,temps,zds,press,obs,watm,satm,mwave,mflux,cdbsloc):
        self.inpath = inpath
        self.outpath = outpath
        self.plotfigs = plotfigs
        self.tags = tags
        self.humids = humids
        self.temps = temps
        self.zds = zds
        self.press = press
        self.obses = obs
        self.watm = watm
        self.satm = satm
        self.mwave0 = mwave
        self.mflux0 = mflux
        self.nights = nights
        self.cdbsloc = cdbsloc


class tagstuffs:

    def __init__(self,night,watm,satm,a0contwave,continuum,ip):
        self.a0contwave = a0contwave
        self.continuum = continuum
        self.watm = watm
        self.satm = satm
        self.night = night
        self.ip = ip
