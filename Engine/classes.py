
# These classes are all for convenient passing of variables between the main code, the telluric spectrum fitter, and the model function.

import numpy as np

class fitobjs:


    def __init__(self,s, x, u,continuum,watm_in,satm_in,mflux_in,mwave_in,mask,masterbeam,CRmask):
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
        self.CRmask = CRmask

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
                                            6: np.array([-0.0000014645, +0.0017169700, +3.3011855894]),
                                            13: np.array([-0.0000012994, +0.0006988905, +3.0840663439]),
                                            14: np.array([-0.0000009176, +0.0018138177, +2.2197275274]),
                                            16: np.array([-0.0000013663, +0.0029909793, +1.3959663664]),
                                            21: np.array([-0.0000013358, +0.0024666477, +1.5425597030]),
                                            22: np.array([-0.0000013951, +0.0029079030, +1.2131742554]),
                                            2: np.array([-0.0000011309, +0.0018503589, +3.2244489850]),
                                            3: np.array([-0.0000006872, +0.0005139971, +3.5055326270]),
                                            4: np.array([-0.0000018218, +0.0026406191, +3.0319937721]),
                                            20: np.array([-0.0000005639, +0.0009729409, +2.0726215568])},
                                        'B':{
                                            6: np.array([-0.0000017946, +0.0022191895, +3.7021645256]),
                                            13: np.array([-0.0000018206, +0.0013364410, +3.3651564498]),
                                            14: np.array([-0.0000009673, +0.0015895234, +2.9497750941]),
                                            16: np.array([-0.0000011332, +0.0020813780, +2.4512905616]),
                                            21: np.array([-0.0000015200, +0.0026951579, +1.8089261830]),
                                            22: np.array([-0.0000017253, +0.0035481206, +1.2235501417]),
                                            2: np.array([-0.0000012003, +0.0017890884, +3.9157171749]),
                                            3: np.array([-0.0000007437, +0.0004339623, +4.1766179383]),
                                            4: np.array([-0.0000021769, +0.0031802110, +3.3217305521]),
                                            20: np.array([-0.0000008535, +0.0014778279, +2.2191974294])}
                                        },
                                    'K':{
                                        'A':{
                                            3: np.array([+0.0000000000, -0.0013882852, +2.6491516427]),
                                            4: np.array([+0.0000005570, -0.0014798964, +3.0322123928]),
                                            5: np.array([+0.0000007815, -0.0018286569, +3.1884587077]),
                                            6: np.array([+0.0000000223, -0.0002438801, +2.5610056734]) },
                                        'B':{
                                            3: np.array([+0.0000000000, -0.0015962816, +3.1179187408]),
                                            4: np.array([+0.0000006349, -0.0021247257, +3.7189580576]),
                                            5: np.array([+0.0000009365, -0.0026350831, +3.9078691486]),
                                            6: np.array([+0.0000003213, -0.0013484977, +3.3906863460]) }
                                        }
                                    }

        self.ips_loosemount_pars = { 'H':{
                                        'A':{
                                            6: np.array([-0.0000016233, +0.0020134037, +3.3487332787]),
                                            13: np.array([-0.0000016702, +0.0011352918, +3.1532836476]),
                                            14: np.array([-0.0000007901, +0.0014705165, +2.5939805593]),
                                            16: np.array([-0.0000014554, +0.0032034428, +1.4193023750]),
                                            21: np.array([-0.0000013562, +0.0024241969, +1.7204152786]),
                                            22: np.array([-0.0000015438, +0.0032647489, +1.1711140400]),
                                            2: np.array([-0.0000009579, +0.0016241564, +3.3749190235]),
                                            3: np.array([-0.0000005287, +0.0002142487, +3.7220784200]),
                                            4: np.array([-0.0000019714, +0.0028921005, +3.0723491330]),
                                            20: np.array([-0.0000009282, +0.0018432727, +1.7490182162])},
                                        'B':{
                                            6: np.array([-0.0000017564, +0.0020515363, +4.0061364082]),
                                            13: np.array([-0.0000020047, +0.0015672027, +3.4707598384]),
                                            14: np.array([-0.0000008056, +0.0012756590, +3.2498824876]),
                                            16: np.array([-0.0000013302, +0.0026174646, +2.2733250373]),
                                            21: np.array([-0.0000015477, +0.0028004919, +1.8687160412]),
                                            22: np.array([-0.0000018236, +0.0037493865, +1.2926912294]),
                                            2: np.array([-0.0000009384, +0.0013175098, +4.1986139959]),
                                            3: np.array([-0.0000006724, +0.0003753155, +4.2393330578]),
                                            4: np.array([-0.0000022722, +0.0033942468, +3.3346821221]),
                                            20: np.array([-0.0000014013, +0.0027124406, +1.7389081440])}
                                        },
                                    'K':{
                                        'A':{
                                            3: np.array([+0.0000000000, -0.0032224396, +5.8539233388]),
                                            4: np.array([+0.0000009369, -0.0042696048, +6.7611052624]),
                                            5: np.array([+0.0000010962, -0.0044543760, +6.7225421397]),
                                            6: np.array([+0.0000003577, -0.0027097128, +5.7789655657]) },
                                        'B':{
                                            3: np.array([+0.0000000000, -0.0033096773, +6.6865750935]),
                                            4: np.array([+0.0000004442, -0.0034871860, +7.1654729260]),
                                            5: np.array([+0.0000006937, -0.0038708805, +7.1443293518]),
                                            6: np.array([-0.0000005376, -0.0011557174, +5.8402964783]) }
                                        }
                                    }


# old H tight~ 0.0149323 , 0.02518648, 0.00716132, 0.00592511, 0.04312596, 0.0075384
        #                                                6           14          16         21          22
        self.methodvariance_tight = { 'H': np.array([0.01145274, 0.00897732, 0.00569098, 0.05562397, 0.0049196]),
        #                                                 3          4           5           6
                                      'K': np.array([0.00184696, 0.0016166, 0.00291516, 0.00139631])
                                    }
        self.methodvariance_loose = { 'H': np.array([np.nan    , np.nan    , np.nan    , np.nan    ]),
                                      'K': np.array([0.0110391, 0.02206625, 0.00439241, 0.00070363])
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
                                    3:  [150, 1350],
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
                                            6: np.array([-0.0000014645, +0.0017169700, +3.3011855894]),
                                            13: np.array([-0.0000012994, +0.0006988905, +3.0840663439]),
                                            14: np.array([-0.0000009176, +0.0018138177, +2.2197275274]),
                                            16: np.array([-0.0000013663, +0.0029909793, +1.3959663664]),
                                            21: np.array([-0.0000013358, +0.0024666477, +1.5425597030]),
                                            22: np.array([-0.0000013951, +0.0029079030, +1.2131742554]),
                                            2: np.array([-0.0000011309, +0.0018503589, +3.2244489850]),
                                            3: np.array([-0.0000006872, +0.0005139971, +3.5055326270]),
                                            4: np.array([-0.0000018218, +0.0026406191, +3.0319937721]),
                                            20: np.array([-0.0000005639, +0.0009729409, +2.0726215568])},
                                        'B':{
                                            6: np.array([-0.0000017946, +0.0022191895, +3.7021645256]),
                                            13: np.array([-0.0000018206, +0.0013364410, +3.3651564498]),
                                            14: np.array([-0.0000009673, +0.0015895234, +2.9497750941]),
                                            16: np.array([-0.0000011332, +0.0020813780, +2.4512905616]),
                                            21: np.array([-0.0000015200, +0.0026951579, +1.8089261830]),
                                            22: np.array([-0.0000017253, +0.0035481206, +1.2235501417]),
                                            2: np.array([-0.0000012003, +0.0017890884, +3.9157171749]),
                                            3: np.array([-0.0000007437, +0.0004339623, +4.1766179383]),
                                            4: np.array([-0.0000021769, +0.0031802110, +3.3217305521]),
                                            20: np.array([-0.0000008535, +0.0014778279, +2.2191974294])}
                                        },
                                    'K':{
                                        'A':{
                                            3: np.array([+0.0000000000, -0.0013882852, +2.6491516427]),
                                            4: np.array([+0.0000005570, -0.0014798964, +3.0322123928]),
                                            5: np.array([+0.0000007815, -0.0018286569, +3.1884587077]),
                                            6: np.array([+0.0000000223, -0.0002438801, +2.5610056734]) },
                                        'B':{
                                            3: np.array([+0.0000000000, -0.0015962816, +3.1179187408]),
                                            4: np.array([+0.0000006349, -0.0021247257, +3.7189580576]),
                                            5: np.array([+0.0000009365, -0.0026350831, +3.9078691486]),
                                            6: np.array([+0.0000003213, -0.0013484977, +3.3906863460]) }
                                        }
                                    }

        self.ips_loosemount_pars = { 'H':{
                                        'A':{
                                            6: np.array([-0.0000016233, +0.0020134037, +3.3487332787]),
                                            13: np.array([-0.0000016702, +0.0011352918, +3.1532836476]),
                                            14: np.array([-0.0000007901, +0.0014705165, +2.5939805593]),
                                            16: np.array([-0.0000014554, +0.0032034428, +1.4193023750]),
                                            21: np.array([-0.0000013562, +0.0024241969, +1.7204152786]),
                                            22: np.array([-0.0000015438, +0.0032647489, +1.1711140400]),
                                            2: np.array([-0.0000009579, +0.0016241564, +3.3749190235]),
                                            3: np.array([-0.0000005287, +0.0002142487, +3.7220784200]),
                                            4: np.array([-0.0000019714, +0.0028921005, +3.0723491330]),
                                            20: np.array([-0.0000009282, +0.0018432727, +1.7490182162])},
                                        'B':{
                                            6: np.array([-0.0000017564, +0.0020515363, +4.0061364082]),
                                            13: np.array([-0.0000020047, +0.0015672027, +3.4707598384]),
                                            14: np.array([-0.0000008056, +0.0012756590, +3.2498824876]),
                                            16: np.array([-0.0000013302, +0.0026174646, +2.2733250373]),
                                            21: np.array([-0.0000015477, +0.0028004919, +1.8687160412]),
                                            22: np.array([-0.0000018236, +0.0037493865, +1.2926912294]),
                                            2: np.array([-0.0000009384, +0.0013175098, +4.1986139959]),
                                            3: np.array([-0.0000006724, +0.0003753155, +4.2393330578]),
                                            4: np.array([-0.0000022722, +0.0033942468, +3.3346821221]),
                                            20: np.array([-0.0000014013, +0.0027124406, +1.7389081440])}
                                        },
                                    'K':{
                                        'A':{
                                            3: np.array([+0.0000000000, -0.0032224396, +5.8539233388]),
                                            4: np.array([+0.0000009369, -0.0042696048, +6.7611052624]),
                                            5: np.array([+0.0000010962, -0.0044543760, +6.7225421397]),
                                            6: np.array([+0.0000003577, -0.0027097128, +5.7789655657]) },
                                        'B':{
                                            3: np.array([+0.0000000000, -0.0033096773, +6.6865750935]),
                                            4: np.array([+0.0000004442, -0.0034871860, +7.1654729260]),
                                            5: np.array([+0.0000006937, -0.0038708805, +7.1443293518]),
                                            6: np.array([-0.0000005376, -0.0011557174, +5.8402964783]) }
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
                                    3:  [150, 1350],
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
