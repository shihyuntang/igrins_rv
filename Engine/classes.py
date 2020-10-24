
# These classes are all for convenient passing of variables between the main code, the telluric spectrum fitter, and the model function.

import numpy as np

class fitobjs:


    def __init__(self,s, x, u,continuum,watm_in,satm_in,mflux_in,mwave_in,mask,masterbeam):
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
                                        'A':{
                                            6: np.array([-0.0000000907, -0.0013918040, +4.9063986084]),
                                            13: np.array([-0.0000000075, -0.0008967200, +3.5674466101]),
                                            14: np.array([-0.0000000234, -0.0007437240, +3.9874010987]),
                                            16: np.array([+0.0000000097, -0.0004680514, +3.4326580933]),
                                            21: np.array([-0.0000001105, -0.0003858732, +3.0295511173]),
                                            22: np.array([+0.0000000127, -0.0003942108, +2.9929330392]) },
                                        'B':{
                                            6: np.array([-0.0000007805, -0.0001433776, +4.9683874066]),
                                            13: np.array([-0.0000008089, -0.0000286379, +3.7677794483]),
                                            14: np.array([-0.0000000122, -0.0010996423, +4.7620526343]),
                                            16: np.array([-0.0000000322, -0.0006017614, +3.9494749224]),
                                            21: np.array([-0.0000002096, -0.0001923039, +3.2665745549]),
                                            22: np.array([-0.0000001010, -0.0002327322, +3.2215225286]) }
                                        },
                                    'K':{
                                        'A':{
                                            3: np.array([-0.0000000179, -0.0008307313, +2.7041105969]),
                                            4: np.array([+0.0000003988, -0.0011221133, +2.8640474133]),
                                            5: np.array([+0.0000001794, -0.0004068603, +2.4287976642]),
                                            6: np.array([-0.0000000056, -0.0001937653, +2.5525097413]) },
                                        'B':{
                                            3: np.array([-0.0000000285, -0.0009487779, +2.9669313124]),
                                            4: np.array([+0.0000002823, -0.0014026499, +3.3962567447]),
                                            5: np.array([+0.0000009472, -0.0026624204, +3.9352693339]),
                                            6: np.array([+0.0000002445, -0.0011871198, +3.3402434648]) }
                                        }
                                    }

        self.ips_loosemount_pars = { 'H':{
                                        'A':{
                                            6: np.array([-0.0000000111, -0.0014319061, +5.0406987618]),
                                            13: np.array([-0.0000004056, +0.0000687306, +3.3113266538]),
                                            14: np.array([-0.0000000593, -0.0005520260, +3.9579877110]),
                                            16: np.array([-0.0000000296, -0.0003538731, +3.5039226324]),
                                            21: np.array([-0.0000000156, -0.0005799910, +3.2513995743]),
                                            22: np.array([-0.0000000274, -0.0003179928, +3.1229337196]) },
                                        'B':{
                                            6: np.array([-0.0000009610, +0.0004310319, +4.7197824011]),
                                            13: np.array([-0.0000009167, +0.0000182371, +3.9320321354]),
                                            14: np.array([+0.0000000098, -0.0010455767, +4.8384723652]),
                                            16: np.array([-0.0000000210, -0.0006543248, +4.1910383559]),
                                            21: np.array([-0.0000000647, -0.0005115389, +3.5571222253]),
                                            22: np.array([-0.0000000465, -0.0004679541, +3.5750046829]) }
                                        },
                                    'K':{
                                        'A':{
                                            3: np.array([-0.0000000247, -0.0023111875, +5.6473632546]),
                                            4: np.array([+0.0000000236, -0.0021640586, +5.6865774873]),
                                            5: np.array([-0.0000000252, -0.0017115857, +5.2320755865]),
                                            6: np.array([-0.0000000517, -0.0018411894, +5.4370535269]) },
                                        'B':{
                                            3: np.array([-0.0000000515, -0.0028328039, +6.8208414627]),
                                            4: np.array([+0.0000000641, -0.0026320115, +6.7329753819]),
                                            5: np.array([+0.0000003191, -0.0028981011, +6.5658524644]),
                                            6: np.array([-0.0000000224, -0.0024526607, +6.5690299081]) }
                                        }
                                    }


# old H tight~ 0.0149323 , 0.02518648, 0.00716132, 0.00592511, 0.04312596, 0.0075384
        #                                                6           13          14          16         21          22
        self.methodvariance_tight = { 'H': np.array([0.02356771, 0.02582711, 0.00694575, 0.00571171, 0.0675675, 0.00698738]),
        #                                                 3          4           5           6
                                      'K': np.array([0.00359593, 0.00159178, 0.00290483, 0.00131022])
                                    }
        self.methodvariance_loose = { 'H': np.array([np.nan    , np.nan    , np.nan    , np.nan    ]),
                                      'K': np.array([0.01796749, 0.01938923, 0.00467536, 0.0012493])
                                    }

        self.bound_cut_dic ={ 'H':{
                                    6: [425, 200],
                                    10: [250, 150],
                                    11: [600, 150],
                                    #13: [200, 600],
                                    13: [250, 750],
                                    14: [750, 100],
                                    16: [530, 100],
                                    17: [1000,100],
                                    20: [500, 150],
                                    21: [200, 150],
                                    22: [200, 150]},
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
        self.maskdict = maskdict
        self.ips_tightmount_pars = { 'H':{
                                        'A':{
                                            6: np.array([-0.0000000907, -0.0013918040, +4.9063986084]),
                                            13: np.array([-0.0000000075, -0.0008967200, +3.5674466101]),
                                            14: np.array([-0.0000000234, -0.0007437240, +3.9874010987]),
                                            16: np.array([+0.0000000097, -0.0004680514, +3.4326580933]),
                                            21: np.array([-0.0000001105, -0.0003858732, +3.0295511173]),
                                            22: np.array([+0.0000000127, -0.0003942108, +2.9929330392]) },
                                        'B':{
                                            6: np.array([-0.0000007805, -0.0001433776, +4.9683874066]),
                                            13: np.array([-0.0000008089, -0.0000286379, +3.7677794483]),
                                            14: np.array([-0.0000000122, -0.0010996423, +4.7620526343]),
                                            16: np.array([-0.0000000322, -0.0006017614, +3.9494749224]),
                                            21: np.array([-0.0000002096, -0.0001923039, +3.2665745549]),
                                            22: np.array([-0.0000001010, -0.0002327322, +3.2215225286]) }
                                        },
                                    'K':{
                                        'A':{
                                            3: np.array([-0.0000000179, -0.0008307313, +2.7041105969]),
                                            4: np.array([+0.0000003988, -0.0011221133, +2.8640474133]),
                                            5: np.array([+0.0000001794, -0.0004068603, +2.4287976642]),
                                            6: np.array([-0.0000000056, -0.0001937653, +2.5525097413]) },
                                        'B':{
                                            3: np.array([-0.0000000285, -0.0009487779, +2.9669313124]),
                                            4: np.array([+0.0000002823, -0.0014026499, +3.3962567447]),
                                            5: np.array([+0.0000009472, -0.0026624204, +3.9352693339]),
                                            6: np.array([+0.0000002445, -0.0011871198, +3.3402434648]) }
                                        }
                                    }

        self.ips_loosemount_pars = { 'H':{
                                        'A':{
                                            6: np.array([-0.0000000111, -0.0014319061, +5.0406987618]),
                                            13: np.array([-0.0000004056, +0.0000687306, +3.3113266538]),
                                            14: np.array([-0.0000000593, -0.0005520260, +3.9579877110]),
                                            16: np.array([-0.0000000296, -0.0003538731, +3.5039226324]),
                                            21: np.array([-0.0000000156, -0.0005799910, +3.2513995743]),
                                            22: np.array([-0.0000000274, -0.0003179928, +3.1229337196]) },
                                        'B':{
                                            6: np.array([-0.0000009610, +0.0004310319, +4.7197824011]),
                                            13: np.array([-0.0000009167, +0.0000182371, +3.9320321354]),
                                            14: np.array([+0.0000000098, -0.0010455767, +4.8384723652]),
                                            16: np.array([-0.0000000210, -0.0006543248, +4.1910383559]),
                                            21: np.array([-0.0000000647, -0.0005115389, +3.5571222253]),
                                            22: np.array([-0.0000000465, -0.0004679541, +3.5750046829]) }
                                        },
                                    'K':{
                                        'A':{
                                            3: np.array([-0.0000000247, -0.0023111875, +5.6473632546]),
                                            4: np.array([+0.0000000236, -0.0021640586, +5.6865774873]),
                                            5: np.array([-0.0000000252, -0.0017115857, +5.2320755865]),
                                            6: np.array([-0.0000000517, -0.0018411894, +5.4370535269]) },
                                        'B':{
                                            3: np.array([-0.0000000515, -0.0028328039, +6.8208414627]),
                                            4: np.array([+0.0000000641, -0.0026320115, +6.7329753819]),
                                            5: np.array([+0.0000003191, -0.0028981011, +6.5658524644]),
                                            6: np.array([-0.0000000224, -0.0024526607, +6.5690299081]) }
                                        }
                                    }


        self.bound_cut_dic ={ 'H':{
                                    6: [425, 200],
                                    10: [250, 150],
                                    11: [600, 150],
                                    #13: [200, 600],
                                    13: [250, 750],
                                    14: [750, 100],
                                    16: [530, 100],
                                    17: [1000,100],
                                    20: [500, 150],
                                    21: [200, 150],
                                    22: [200, 150]},
                               'K':{
                                    13: [200, 400],
                                    14: [200, 400]}
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
