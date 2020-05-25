
# These classes are all for convenient passing of variables between the main code, the telluric spectrum fitter, and the model function.

import numpy as np

class fitobjs:

    def __init__(self,s, x, u,watm_in,satm_in,uatm_in,mflux_in,mwave_in):
        self.s = s
        self.x = x
        self.u = u
        self.uatm_in = uatm_in
        self.watm_in = watm_in
        self.satm_in = satm_in
        self.mflux_in = mflux_in
        self.mwave_in = mwave_in

class inparams:

    def __init__(self,inpath,outpath,initvsini,vsinivary,plotfigs,initguesses,bvcs,tagsA,tagsB,nights,mwave,mflux,a0dict,xbounddict):
        self.inpath      = inpath
        self.outpath     = outpath
        self.initvsini   = initvsini
        self.vsinivary   = vsinivary
        self.plotfigs    = plotfigs
        self.initguesses = initguesses
        self.bvcs        = bvcs
        self.tagsA       = tagsA
        self.tagsB       = tagsB
        self.mwave0      = mwave
        self.mflux0      = mflux
        self.nights      = nights
        self.a0dict      = a0dict
        self.xbounddict  = xbounddict
        self.ips_tightmount_pars = { 'H':{
                                             2: np.array([-0.00000083, +0.00141598, 3.6143533]),
                                             3: np.array([-0.00000042, -0.00009824, 3.9964264]),#CITau + GJ281 full
                                             4: np.array([-0.00000164, +0.00226626, 3.3720807]),#CITau + GJ281 full
                                             5: np.array([+0.00000000, -0.00061770, 4.0095721]),#CITau + GJ281 full
                                             6: np.array([-0.00000035, -0.00074569, 4.8377666]),#CITau + GJ281 full
                                            10: np.array([-0.00000051, +0.00007716, 4.1389145]),#CITau + GJ281 full
                                            11: np.array([-0.00000000, -0.00045612, 3.6545044]),
                                            13: np.array([-0.00000055, -0.00057129, 3.8362119]),#CITau + GJ281 full
                                            14: np.array([-0.00000001, -0.00087136, 4.2944599]),#CITau + GJ281 full
                                            16: np.array([-0.00000000, -0.00052450, 3.6611907]),#CITau + GJ281 full
                                            17: np.array([-0.00000000, -0.00000074, 2.6779583]),
                                            20: np.array([+0.00000001, -0.00041767, 3.0446818]),
                                            21: np.array([+0.00000001, -0.00059766, 3.3326337]),#CITau + GJ281 full
                                            22: np.array([+0.00000000, -0.00036886, 3.2415346]) },
                                    'K':{    2: np.array([-0.00000019, +0.00028060, 1.3383606]),#CITau + GJ281 full
                                             3: np.array([-0.00000014, -0.00066076, 2.9978030]),#CITau + GJ281 full
                                             4: np.array([+0.00000014, -0.00098005, 3.3793699]),#CITau + GJ281 full
                                             5: np.array([+0.00000032, -0.00096880, 3.1177634]),#CITau + GJ281 full
                                             6: np.array([+0.00000006, -0.00066340, 3.2028759]),#CITau + GJ281 full
                                             7: np.array([+0.00000013, -0.00069198, 3.4463934]),
                                             8: np.array([+0.00000042, -0.00102179, 3.037313]),#CITau + GJ281 full
                                            10: np.array([+0.00000001, +0.00003481, 2.1329333]),#CITau + GJ281 full
                                            11: np.array([-0.00000021, +0.00020740, 2.5122968]),#CITau + GJ281 full
                                            12: np.array([+0.00000000, -0.00070081, 2.2782001]),
                                            13: np.array([+0.00000000, -0.00049637, 1.8553702]),
                                            14: np.array([+0.00000040, -0.00104277, 2.8110920]),
                                            16: np.array([-0.00000000, +0.00063948, 1.7424173]) }
                                    }
        self.ips_loosemount_pars = { 'H':{
                                             2:  np.array([-0.00000098, +0.00179662, 3.6487848]),
                                             3: np.array([-0.00000040, -0.00020271, 4.3306463]),#CITau + GJ281 full
                                             4: np.array([-0.00000171, +0.00225548, 3.7094981]),#CITau + GJ281 full
                                             5: np.array([-0.00000012, -0.00039063, 4.1635584]),#CITau + GJ281 full
                                             6: np.array([-0.00000040, -0.00058547, 4.8982556]),#CITau + GJ281 full
                                            10: np.array([-0.00000062, +0.00025086, 4.3092912]),#CITau + GJ281 full
                                            11: np.array([-0.00000001, -0.00043510, 4.0044980]),
                                            13: np.array([-0.00000065, -0.00049497, 4.0052192]),#CITau + GJ281 full
                                            14: np.array([-0.00000026, -0.00038912, 4.3156883]),#CITau + GJ281 full
                                            16: np.array([-0.00000000, -0.00062674, 3.9625779]),#CITau + GJ281 full
                                            17: np.array([-0.00000000, -0.00054249, 3.4089812]),
                                            20: np.array([+0.00000000, -0.00050703, 3.3124597]),
                                            21: np.array([+0.00000000, -0.00059269, 3.4646404]),#CITau + GJ281 full
                                            22: np.array([-0.00000000, -0.00037943, 3.4248205]) },
                                    'K':{    2: np.array([-0.00000030, -0.00136354, 5.2787831]),#CITau + GJ281 full
                                             3: np.array([-0.00000064, -0.00055100, 5.0549681]),#CITau + GJ281 full
                                             4: np.array([-0.00000029, -0.00130082, 5.5659733]),#CITau + GJ281 full
                                             5: np.array([-0.00000051, -0.00034471, 4.7191411]),#CITau + GJ281 full
                                             6: np.array([-0.00000116, +0.00073074, 4.4859013]),#CITau + GJ281 full
                                             7: np.array([-0.00000027, -0.00040471, 4.2938356]),
                                             8: np.array([-0.00000022, -0.00072976, 4.4175125]),#CITau + GJ281 full
                                            10: np.array([-0.00000023, -0.00107077, 4.3597473]),#CITau + GJ281 full
                                            11: np.array([-0.00000007, -0.00122038, 4.431836]),#CITau + GJ281 full
                                            12: np.array([+0.00000006, -0.00120924, 3.2901135]),
                                            13: np.array([-0.00000000, -0.00116862, 3.0264358]),
                                            14: np.array([+0.00000017, -0.00142531, 4.0935727]),
                                            16: np.array([-0.00000000, -0.00008538, 2.6588594]) }
                                    }
        #
        # self.ips_tightmount = {2: np.array([2.04037987, 2.44724324, 1.94710219, 2.25831301, 2.04870909,
        #                            2.54918255]), 3: np.array([2.5835808 , 2.51815444, 2.54622661, 2.33003733, 2.21602917,
        #                            1.6183255 ]), 4: np.array([2.60463925, 2.50289512, 2.18010336, 2.08725543, 2.0447105 ,
        #                            2.68329089]), 5: np.array([2.62590674, 2.34995633, 2.26398339, 2.10716032, 2.16699697,
        #                            2.60612536]), 6: np.array([2.58312563, 2.4476741 , 2.28593311, 2.18518033, 2.50510778,
        #                            2.63379758])}
        # self.ips_loosemount = {2: np.array([5.23186917, 5.24900976, 3.96609885, 3.85262363, 3.11525618,
        #                            3.6394157 ]), 3: np.array([5.68122268, 4.94600199, 4.42893859, 3.97849734, 3.1649909 ,
        #                            2.35128276]), 4: np.array([5.49237893, 5.02570037, 3.91499657, 3.44131085, 3.04847074,
        #                            3.3096322 ]), 5: np.array([5.63882553, 4.60165159, 4.01001552, 3.40237606, 3.06376238,
        #                            2.98584427]), 6: np.array([5.25760553, 4.56458826, 3.93430453, 3.30376016, 3.11607041,
        #                            3.0286326 ])}
        # self.ips_tightmountFunc = {2: np.array([ 2.25378271e-09, -6.44641469e-06,  5.56041474e-03,  7.72955295e-01]),
        #                             3: np.array([-1.03836176e-09,  2.28983215e-06, -1.75341703e-03,  2.97931406e+00]),
        #                             4: np.array([ 2.20182527e-09, -5.41154039e-06,  3.32290632e-03,  2.00177689e+00]),
        #                             5: np.array([ 1.00979196e-09, -1.94770701e-06,  3.99506340e-04,  2.68272193e+00]),
        #                             6: np.array([ 1.52820030e-10,  4.24369854e-07, -1.35635083e-03,  3.04970584e+00])}
        #
        # self.ips_loosemountFunc = {2: np.array([ 4.05939305e-09, -1.11231891e-05,  7.06911774e-03,  3.98215719e+00]),
        #                             3: np.array([-1.38630177e-09,  3.80492542e-06, -5.48547631e-03,  7.28241681e+00]),
        #                             4: np.array([ 2.62862802e-09, -6.25435468e-06,  1.75570160e-03,  5.62385973e+00]),
        #                             5: np.array([-7.39811208e-11,  1.81851999e-06, -5.51876063e-03,  7.46185799e+00]),
        #                             6: np.array([ 8.12015898e-10, -1.16758952e-06, -2.22011098e-03,  6.23343180e+00])}
        #
        #
        # self.chunkweights_loosemount = {2: np.array([0.13992335, 0.00548394, 0.37635348, 0.21547168, 0.03589399, 0.22687356]),
        #                                 3: np.array([0.20653672, 0.16309499, 0.06580472, 0.04768457, 0.31613504, 0.20074395]),
        #                                 4: np.array([0.06278375, 0.18963667, 0.05826427, 0.25202664, 0.27911588, 0.1581728 ]),
        #                                 5: np.array([0.05326762, 0.02296659, 0.03711518, 0.33631157, 0.3878665 , 0.16247254]),
        #                                 6: np.array([0.03406971, 0.19150636, 0.15882698, 0.09657096, 0.28659474, 0.23243126])}
        #
        # self.chunkweights_tightmount = {2: np.array([0.21525525, 0.01333442, 0.61514991, 0.05957497, 0.01979944, 0.07688601]),
        #                                 3: np.array([0.21843997, 0.26091404, 0.13681933, 0.06954126, 0.13812044, 0.17616496]),
        #                                 4: np.array([0.1526929 , 0.29268455, 0.11404597, 0.19192788, 0.19955444, 0.04909426]),
        #                                 5: np.array([0.13061788, 0.07999278, 0.055955  , 0.3347575 , 0.22324851, 0.17542833]),
        #                                 6: np.array([0.05020937, 0.30509215, 0.33481064, 0.08170708, 0.16561447, 0.06256628])}
        self.methodvariance_tight = np.array([0.00282905, 0.00090454, 0.00124326, 0.00210822, 0.00140632])
        self.methodvariance_loose = np.array([0.00883368, 0.00840557, 0.00534123, 0.00385086, 0.00219453])
        #






class inparamsA0:

    def __init__(self,inpath,outpath,plotfigs,tags,nights,humids,temps,zds,press,obs,watm,satm,mwave,mflux,cdbsloc,xbounddict):
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
        self.xbounddict  = xbounddict
        self.ips_tightmount_pars = { 'H':{
                                             2: np.array([-0.00000083, +0.00141598, 3.6143533]),
                                             3: np.array([-0.00000042, -0.00009824, 3.9964264]),#CITau + GJ281 full
                                             4: np.array([-0.00000164, +0.00226626, 3.3720807]),#CITau + GJ281 full
                                             5: np.array([+0.00000000, -0.00061770, 4.0095721]),#CITau + GJ281 full
                                             6: np.array([-0.00000035, -0.00074569, 4.8377666]),#CITau + GJ281 full
                                            10: np.array([-0.00000051, +0.00007716, 4.1389145]),#CITau + GJ281 full
                                            11: np.array([-0.00000000, -0.00045612, 3.6545044]),
                                            13: np.array([-0.00000055, -0.00057129, 3.8362119]),#CITau + GJ281 full
                                            14: np.array([-0.00000001, -0.00087136, 4.2944599]),#CITau + GJ281 full
                                            16: np.array([-0.00000000, -0.00052450, 3.6611907]),#CITau + GJ281 full
                                            17: np.array([-0.00000000, -0.00000074, 2.6779583]),
                                            20: np.array([+0.00000001, -0.00041767, 3.0446818]),
                                            21: np.array([+0.00000001, -0.00059766, 3.3326337]),#CITau + GJ281 full
                                            22: np.array([+0.00000000, -0.00036886, 3.2415346]) },
                                    'K':{    2: np.array([-0.00000019, +0.00028060, 1.3383606]),#CITau + GJ281 full
                                             3: np.array([-0.00000014, -0.00066076, 2.9978030]),#CITau + GJ281 full
                                             4: np.array([+0.00000014, -0.00098005, 3.3793699]),#CITau + GJ281 full
                                             5: np.array([+0.00000032, -0.00096880, 3.1177634]),#CITau + GJ281 full
                                             6: np.array([+0.00000006, -0.00066340, 3.2028759]),#CITau + GJ281 full
                                             7: np.array([+0.00000013, -0.00069198, 3.4463934]),
                                             8: np.array([+0.00000042, -0.00102179, 3.037313]),#CITau + GJ281 full
                                            10: np.array([+0.00000001, +0.00003481, 2.1329333]),#CITau + GJ281 full
                                            11: np.array([-0.00000021, +0.00020740, 2.5122968]),#CITau + GJ281 full
                                            12: np.array([+0.00000000, -0.00070081, 2.2782001]),
                                            13: np.array([+0.00000000, -0.00049637, 1.8553702]),
                                            14: np.array([+0.00000040, -0.00104277, 2.8110920]),
                                            16: np.array([-0.00000000, +0.00063948, 1.7424173]) }
                                    }
        self.ips_loosemount_pars = { 'H':{
                                             2:  np.array([-0.00000098, +0.00179662, 3.6487848]),
                                             3: np.array([-0.00000040, -0.00020271, 4.3306463]),#CITau + GJ281 full
                                             4: np.array([-0.00000171, +0.00225548, 3.7094981]),#CITau + GJ281 full
                                             5: np.array([-0.00000012, -0.00039063, 4.1635584]),#CITau + GJ281 full
                                             6: np.array([-0.00000040, -0.00058547, 4.8982556]),#CITau + GJ281 full
                                            10: np.array([-0.00000062, +0.00025086, 4.3092912]),#CITau + GJ281 full
                                            11: np.array([-0.00000001, -0.00043510, 4.0044980]),
                                            13: np.array([-0.00000065, -0.00049497, 4.0052192]),#CITau + GJ281 full
                                            14: np.array([-0.00000026, -0.00038912, 4.3156883]),#CITau + GJ281 full
                                            16: np.array([-0.00000000, -0.00062674, 3.9625779]),#CITau + GJ281 full
                                            17: np.array([-0.00000000, -0.00054249, 3.4089812]),
                                            20: np.array([+0.00000000, -0.00050703, 3.3124597]),
                                            21: np.array([+0.00000000, -0.00059269, 3.4646404]),#CITau + GJ281 full
                                            22: np.array([-0.00000000, -0.00037943, 3.4248205]) },
                                    'K':{    2: np.array([-0.00000030, -0.00136354, 5.2787831]),#CITau + GJ281 full
                                             3: np.array([-0.00000064, -0.00055100, 5.0549681]),#CITau + GJ281 full
                                             4: np.array([-0.00000029, -0.00130082, 5.5659733]),#CITau + GJ281 full
                                             5: np.array([-0.00000051, -0.00034471, 4.7191411]),#CITau + GJ281 full
                                             6: np.array([-0.00000116, +0.00073074, 4.4859013]),#CITau + GJ281 full
                                             7: np.array([-0.00000027, -0.00040471, 4.2938356]),
                                             8: np.array([-0.00000022, -0.00072976, 4.4175125]),#CITau + GJ281 full
                                            10: np.array([-0.00000023, -0.00107077, 4.3597473]),#CITau + GJ281 full
                                            11: np.array([-0.00000007, -0.00122038, 4.431836]),#CITau + GJ281 full
                                            12: np.array([+0.00000006, -0.00120924, 3.2901135]),
                                            13: np.array([-0.00000000, -0.00116862, 3.0264358]),
                                            14: np.array([+0.00000017, -0.00142531, 4.0935727]),
                                            16: np.array([-0.00000000, -0.00008538, 2.6588594]) }
                                    }


class orderdict_cla:
    def __init__(self,):
        self.orderdict = { 'H':{1: [1.79350, 1.81560],
                                2: [1.77602, 1.79791],
                                3: [1.75889, 1.78058],
                                4: [1.74211, 1.76360],
                                5: [1.72565, 1.74694],
                                6: [1.70952, 1.73061],
                                7: [1.69811, 1.71032],
                                8: [1.68255, 1.69465],
                                9: [1.66729, 1.67928],
                                10: [1.64928, 1.66785],
                                11: [1.63880, 1.65301],
                                12: [1.62320, 1.63487],
                                13: [1.60548, 1.61957],
                                14: [1.59742, 1.61061],
                                15: [1.58149, 1.59284],
                                16: [1.56694, 1.58328],
                                17: [1.56034, 1.56997],
                                18: [1.54197, 1.55303],
                                19: [1.52926, 1.54022],
                                20: [1.51677, 1.53096],
                                21: [1.50060, 1.51902],
                                22: [1.48857, 1.50682],
                                23: [1.47673, 1.49483],
                                24: [1.46510, 1.48303],
                                25: [1.45366, 1.47143]},
                           'K':{1: [2.45317, 2.48280],
                                2: [2.41999, 2.44928],
                                3: [2.38771, 2.41667],
                                4: [2.35630, 2.38493],
                                5: [2.32573, 2.35402],
                                6: [2.29596, 2.32393],
                                7: [2.26696, 2.29461],
                                8: [2.23870, 2.26603],
                                9: [2.21115, 2.23817],
                                10: [2.18429, 2.21101],
                                11: [2.15891, 2.18451],
                                12: [2.14410, 2.15604],
                                13: [2.10840, 2.12953],
                                14: [2.08326, 2.10624],
                                15: [2.05948, 2.08471],
                                16: [2.03627, 2.06121],
                                17: [2.01358, 2.03825],
                                18: [1.99142, 2.01580],
                                19: [1.96975, 1.99386],
                                20: [1.94857, 1.97241],
                                21: [1.92785, 1.95143],
                                22: [1.90759, 1.93090],
                                23: [1.88777, 1.91082],
                                24: [1.86837, 1.89116],
                                25: [1.84939, 1.87192]}
                            }

class tagstuffs:

    def __init__(self,night,watm,satm,a0contwave,continuum,ip):
        self.a0contwave = a0contwave
        self.continuum = continuum
        self.watm = watm
        self.satm = satm
        self.night = night
        self.ip = ip
