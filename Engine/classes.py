
# These classes are all for convenient passing of variables between the main code, the telluric spectrum fitter, and the model function.

import numpy as np

class fitobjs:


    def __init__(self,s, x, u,continuum,watm_in,satm_in,mflux_in,mwave_in,mask,masterbeam=None):
        self.s = s
        self.x = x
        self.u = u
        self.continuum = continuum
        self.watm_in = watm_in
        self.satm_in = satm_in
        self.mflux_in = mflux_in
        self.mwave_in = mwave_in
        self.mask = mask
        self.masterbeam = masterbeam

class inparams:

    def __init__(self,inpath,outpath,initvsini,vsinivary,plotfigs,initguesses,bvcs,tagsA,tagsB,nights,mwave,mflux,a0dict,xbounddict,maskdict):
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
        self.maskdict = maskdict
        self.ips_tightmount_pars = { 'H':{
                                            2: np.array([-0.00000127, +0.00210607, 3.4481905]),
                                            3: np.array([-0.00000086, +0.00113201, 3.4055482]),
                                            4: np.array([-0.00000173, +0.00253586, 3.3329716]),
                                            6: np.array([-0.00000139, +0.00183241, 3.466488]),
                                            13: np.array([-0.00000053, -0.00002286, 3.5243779]),
                                            14: np.array([-0.00000023, +0.00003049, 3.5175536]),
                                            16: np.array([-0.00000010, -0.00031187, 3.6099983]),
                                            20: np.array([+0.00000005, -0.00057104, 3.1805937]),
                                            21: np.array([-0.00000004, -0.00049309, 3.3365371]),
                                            22: np.array([-0.00000005, -0.00031042, 3.2413131])},
                                    'K':{   2: np.array([-0.00000043, -0.00009260, 2.8901846]),
                                            3: np.array([-0.00000029, -0.00045505, 3.110618]),
                                            4: np.array([-0.00000022, -0.00023683, 3.0734207]),
                                            5: np.array([+0.00000020, -0.00084550, 3.1787218]),
                                            6: np.array([-0.00000007, -0.00049691, 3.2275649]),
                                            7: np.array([+0.00000019, -0.00081209, 3.1828768]),
                                            8: np.array([+0.00000037, -0.00108785, 3.2596467]),
                                            10: np.array([+0.00000003, -0.00065472, 2.926037]),
                                            14: np.array([+0.00000036, -0.00153806, 3.1654484]),
                                            16: np.array([+0.00000079, -0.00127711, 2.8345052]) }
                                    }

        self.ips_loosemount_pars = { 'H':{
                                            2: np.array([-0.00000097, +0.00180228, 3.4908878]),
                                            3: np.array([-0.00000083, +0.00098296, 3.6176389]),
                                            4: np.array([-0.00000179, +0.00256885, 3.422362]),
                                            6: np.array([-0.00000137, +0.00197794, 3.3224567]),
                                            13: np.array([-0.00000044, -0.00007479, 3.6248183]),
                                            14: np.array([-0.00000003, -0.00027106, 3.6804099]),
                                            16: np.array([-0.00000001, -0.00036968, 3.6468611]),
                                            20: np.array([+0.00000000, -0.00051203, 3.277703]),
                                            21: np.array([+0.00000000, -0.00053837, 3.3633831]),
                                            22: np.array([-0.00000007, -0.00018335, 3.1698995]) },
                                    'K':{   2: np.array([-0.00000100, +0.00133316, 3.4675222]),
                                            3: np.array([-0.00000068, +0.00070342, 3.7342004]),
                                            4: np.array([-0.00000093, +0.00126044, 3.6353763]),
                                            5: np.array([-0.00000049, +0.00036114, 3.6458918]),
                                            6: np.array([-0.00000044, +0.00024862, 3.7530547]),
                                            7: np.array([-0.00000046, +0.00033508, 3.7023568]),
                                            8: np.array([-0.00000008, -0.00049619, 3.7524707]),
                                            10: np.array([+0.00000000, -0.00077256, 3.6673045]),
                                            14: np.array([+0.00000000, -0.00104468, 3.760597]),
                                            16: np.array([+0.00000029, -0.00070009, 2.9448922]) }
                                    }


# old H tight~ 0.0149323 , 0.02518648, 0.00716132, 0.00592511, 0.04312596, 0.0075384 
        #                                                6           13          14          16         21          22
        self.methodvariance_tight = { 'H': np.array([0.02356771, 0.02582711, 0.00694575, 0.00571171, 0.0675675, 0.00698738]),
        #                                                 3          4           5           6
                                      'K': np.array([0.00403419, 0.00227922, 0.0035236 , 0.00150486])
                                    }
        self.methodvariance_loose = { 'H': np.array([np.nan    , np.nan    , np.nan    , np.nan    ]),
                                      'K': np.array([0.01980407, 0.02014976, 0.00789618, 0.00024278])
                                    }

        self.bound_cut_dic = { 'H':{
                                    10: [250, 150],
                                    11: [600, 150],
                                    13: [200, 600],
                                    14: [700, 100],
                                    16: [400, 100],
                                    17: [1000,100],
                                    20: [500, 150]},
                               'K':{
                                    13: [200, 400],
                                    14: [200, 400]}
                            }

class inparamsA0:

    def __init__(self,inpath,outpath,plotfigs,tags,nights,humids,temps,zds,press,obs,watm,satm,mwave,mflux,cdbsloc,xbounddict,maskdict):
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
        self.maskdict = maskdict
        self.ips_tightmount_pars = { 'H':{
                                            2: np.array([-0.00000127, +0.00210607, 3.4481905]),
                                            3: np.array([-0.00000086, +0.00113201, 3.4055482]),
                                            4: np.array([-0.00000173, +0.00253586, 3.3329716]),
                                            6: np.array([-0.00000139, +0.00183241, 3.466488]),
                                            13: np.array([-0.00000053, -0.00002286, 3.5243779]),
                                            14: np.array([-0.00000023, +0.00003049, 3.5175536]),
                                            16: np.array([-0.00000010, -0.00031187, 3.6099983]),
                                            20: np.array([+0.00000005, -0.00057104, 3.1805937]),
                                            21: np.array([-0.00000004, -0.00049309, 3.3365371]),
                                            22: np.array([-0.00000005, -0.00031042, 3.2413131])},
                                    'K':{   2: np.array([-0.00000043, -0.00009260, 2.8901846]),
                                            3: np.array([-0.00000029, -0.00045505, 3.110618]),
                                            4: np.array([-0.00000022, -0.00023683, 3.0734207]),
                                            5: np.array([+0.00000020, -0.00084550, 3.1787218]),
                                            6: np.array([-0.00000007, -0.00049691, 3.2275649]),
                                            7: np.array([+0.00000019, -0.00081209, 3.1828768]),
                                            8: np.array([+0.00000037, -0.00108785, 3.2596467]),
                                            10: np.array([+0.00000003, -0.00065472, 2.926037]),
                                            14: np.array([+0.00000036, -0.00153806, 3.1654484]),
                                            16: np.array([+0.00000079, -0.00127711, 2.8345052]) }
                                    }

        self.ips_loosemount_pars = { 'H':{
                                            2: np.array([-0.00000097, +0.00180228, 3.4908878]),
                                            3: np.array([-0.00000083, +0.00098296, 3.6176389]),
                                            4: np.array([-0.00000179, +0.00256885, 3.422362]),
                                            6: np.array([-0.00000137, +0.00197794, 3.3224567]),
                                            13: np.array([-0.00000044, -0.00007479, 3.6248183]),
                                            14: np.array([-0.00000003, -0.00027106, 3.6804099]),
                                            16: np.array([-0.00000001, -0.00036968, 3.6468611]),
                                            20: np.array([+0.00000000, -0.00051203, 3.277703]),
                                            21: np.array([+0.00000000, -0.00053837, 3.3633831]),
                                            22: np.array([-0.00000007, -0.00018335, 3.1698995]) },
                                    'K':{   2: np.array([-0.00000100, +0.00133316, 3.4675222]),
                                            3: np.array([-0.00000068, +0.00070342, 3.7342004]),
                                            4: np.array([-0.00000093, +0.00126044, 3.6353763]),
                                            5: np.array([-0.00000049, +0.00036114, 3.6458918]),
                                            6: np.array([-0.00000044, +0.00024862, 3.7530547]),
                                            7: np.array([-0.00000046, +0.00033508, 3.7023568]),
                                            8: np.array([-0.00000008, -0.00049619, 3.7524707]),
                                            10: np.array([+0.00000000, -0.00077256, 3.6673045]),
                                            14: np.array([+0.00000000, -0.00104468, 3.760597]),
                                            16: np.array([+0.00000029, -0.00070009, 2.9448922]) }
                                    }


        self.bound_cut_dic = { 'H':{
                                    13: [200, 600],
                                    14: [700, 100],
                                    16: [400, 100],
                                    20: [500, 150]},
                               'K':{
                                    14: [150, 300]}
                               }


class orderdict_cla:
    def __init__(self,):
        self.orderdict = { 'H':{1: [1.79350, 1.81560],
                                2: [1.77602, 1.79791],
                                3: [1.75889, 1.78058],
                                4: [1.74211, 1.76360],
                                5: [1.72565, 1.74694],
                                6: [1.70952, 1.73061],
                                # 7: [1.69811, 1.71032],
                                # 8: [1.68255, 1.69465],
                                # 9: [1.66729, 1.67928],
                                10: [1.64928, 1.66785],
                                11: [1.63880, 1.65301],
                                # 12: [1.62320, 1.63487],
                                13: [1.60548, 1.61957],
                                14: [1.59742, 1.61061],
                                # 15: [1.58149, 1.59284],
                                16: [1.56694, 1.58328],
                                17: [1.56034, 1.56997],
                                # 18: [1.54197, 1.55303],
                                # 19: [1.52926, 1.54022],
                                20: [1.51677, 1.53096],
                                21: [1.50060, 1.51902],
                                22: [1.48857, 1.50682],
                                # 23: [1.47673, 1.49483],
                                # 24: [1.46510, 1.48303],
                                25: [1.45366, 1.47143]},
                           'K':{1: [2.45317, 2.48280],
                                2: [2.41999, 2.44928],
                                3: [2.38771, 2.41667],
                                4: [2.35630, 2.38493],
                                5: [2.32573, 2.35402],
                                6: [2.29596, 2.32393],
                                7: [2.26696, 2.29461],
                                8: [2.23870, 2.26603],
                                # 9: [2.21115, 2.23817],
                                10: [2.18429, 2.21101],
                                11: [2.15891, 2.18451],
                                12: [2.14410, 2.15604],
                                13: [2.10840, 2.12953],
                                14: [2.08326, 2.10624],
                                # 15: [2.05948, 2.08471],
                                16: [2.03627, 2.06121],
                                # 17: [2.01358, 2.03825],
                                # 18: [1.99142, 2.01580],
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
