# These classes are all for convenient passing of variables between the main 
# code, the telluric spectrum fitter, and the model function.

import numpy as np

ips_tight = { 
    'H':{
        'A':{
            6: np.array([-0.0000014645, +0.0017169700, +3.3011855894]),
            13: np.array([-0.0000012994, +0.0006988905, +3.0840663439]),
            14: np.array([-0.0000009176, +0.0018138177, +2.2197275274]),
            16: np.array([-0.0000013663, +0.0029909793, +1.3959663664]),
            21: np.array([-0.0000013358, +0.0024666477, +1.5425597030]),
            22: np.array([-0.0000013951, +0.0029079030, +1.2131742554]),
            2: np.array([-0.0000011308, +0.0019048124, +3.1673497542]),
            3: np.array([-0.0000006216, +0.0004114534, +3.4773458923]),
            4: np.array([-0.0000018462, +0.0027342132, +2.9057196779]),
            20: np.array([-0.0000005747, +0.0010462689, +2.0063265625])},
        'B':{
            6: np.array([-0.0000017946, +0.0022191895, +3.7021645256]),
            13: np.array([-0.0000018206, +0.0013364410, +3.3651564498]),
            14: np.array([-0.0000009673, +0.0015895234, +2.9497750941]),
            16: np.array([-0.0000011332, +0.0020813780, +2.4512905616]),
            21: np.array([-0.0000015200, +0.0026951579, +1.8089261830]),
            22: np.array([-0.0000017253, +0.0035481206, +1.2235501417]),
            2: np.array([-0.0000011552, +0.0017634890, +3.8604002255]),
            3: np.array([-0.0000006966, +0.0003965295, +4.0969467142]),
            4: np.array([-0.0000020779, +0.0029574695, +3.3703599880]),
            20: np.array([-0.0000010270, +0.0019515216, +1.9071269518])}
        },
    'K':{
        'A':{
            3: np.array([+0.0000000000, -0.0002572599, +2.5524172079]),
            4: np.array([+0.0000006832, -0.0016426204, +3.1435869145]),
            5: np.array([+0.0000008324, -0.0018784941, +3.2860911563]),
            6: np.array([+0.0000002202, -0.0006810739, +2.7989979812]),
            2: np.array([+0.0000018993, -0.0041545122, +3.2722615480]),
            7: np.array([+0.0000003846, -0.0008248737, +2.7054577962]),
            8: np.array([+0.0000006032, -0.0011831188, +2.7613774604]),
            10: np.array([-0.0000002109, +0.0008841434, +1.2082829518]),
            14: np.array([+0.0000007984, -0.0020449382, +3.1101522478]),
            16: np.array([+0.0000001490, +0.0009528175, +1.1385210594]) },
        'B':{
            3: np.array([+0.0000000000, -0.0007619935, +3.1856884900]),
            4: np.array([+0.0000007733, -0.0022856887, +3.7896224682]),
            5: np.array([+0.0000010009, -0.0026751285, +3.9697999424]),
            6: np.array([+0.0000003561, -0.0014029046, +3.4360799984]),
            2: np.array([+0.0000017922, -0.0041528862, +3.4063864564]),
            7: np.array([+0.0000005729, -0.0016962859, +3.4733681135]),
            8: np.array([+0.0000007786, -0.0020117585, +3.4821592021]),
            10: np.array([-0.0000001630, +0.0003023356, +1.8355214070]),
            14: np.array([+0.0000006373, -0.0022195328, +3.4337657140]),
            16: np.array([+0.0000000959, +0.0007723205, +1.3240275178]) }
        }
    }

ips_loose = { 
    'H':{
        'A':{
            6: np.array([-0.0000016233, +0.0020134037, +3.3487332787]),
            13: np.array([-0.0000016702, +0.0011352918, +3.1532836476]),
            14: np.array([-0.0000007901, +0.0014705165, +2.5939805593]),
            16: np.array([-0.0000014554, +0.0032034428, +1.4193023750]),
            21: np.array([-0.0000013562, +0.0024241969, +1.7204152786]),
            22: np.array([-0.0000015438, +0.0032647489, +1.1711140400]),
            2: np.array([-0.0000007673, +0.0011908119, +3.5574747235]),
            3: np.array([-0.0000004539, +0.0000711467, +3.7569496247]),
            4: np.array([-0.0000017014, +0.0022235883, +3.3469057729]),
            20: np.array([-0.0000011548, +0.0024159347, +1.4899767054])},
        'B':{
            6: np.array([-0.0000017564, +0.0020515363, +4.0061364082]),
            13: np.array([-0.0000020047, +0.0015672027, +3.4707598384]),
            14: np.array([-0.0000008056, +0.0012756590, +3.2498824876]),
            16: np.array([-0.0000013302, +0.0026174646, +2.2733250373]),
            21: np.array([-0.0000015477, +0.0028004919, +1.8687160412]),
            22: np.array([-0.0000018236, +0.0037493865, +1.2926912294]),
            2: np.array([-0.0000009144, +0.0012463162, +4.2410289004]),
            3: np.array([-0.0000005587, +0.0001065024, +4.3728963431]),
            4: np.array([-0.0000021505, +0.0029466667, +3.6128512030]),
            20: np.array([-0.0000014158, +0.0027485459, +1.7502302932])}
        },
    'K':{
        'A':{
            3: np.array([+0.0000000000, -0.0024074076, +6.2739343639]),
            4: np.array([+0.0000011138, -0.0045865523, +7.0065081737]),
            5: np.array([+0.0000012423, -0.0046623467, +6.8935184675]),
            6: np.array([+0.0000006996, -0.0035210592, +6.2893198778]),
            2: np.array([+0.0000000173, -0.0019826425, +4.8330715369]),
            7: np.array([+0.0000008492, -0.0035313488, +5.9016867272]),
            8: np.array([+0.0000010973, -0.0038178120, +5.6832992738]),
            10: np.array([+0.0000002443, -0.0017360668, +4.1415859150]),
            14: np.array([+0.0000015291, -0.0039022985, +4.6823412563]),
            16: np.array([+0.0000003064, -0.0001454813, +2.0998323933]) },
        'B':{
            3: np.array([+0.0000000000, -0.0030490721, +7.4604027469]),
            4: np.array([+0.0000005150, -0.0035314345, +7.2643757724]),
            5: np.array([+0.0000007942, -0.0039583272, +7.2580505977]),
            6: np.array([+0.0000004405, -0.0035119336, +7.2082659894]),
            2: np.array([-0.0000001693, -0.0021423525, +5.9778102108]),
            7: np.array([+0.0000006778, -0.0035919133, +6.7264063526]),
            8: np.array([+0.0000010518, -0.0041659367, +6.6036171357]),
            10: np.array([+0.0000000258, -0.0016619486, +4.7519124026]),
            14: np.array([-0.0000007761, -0.0001818842, +3.9448840447]),
            16: np.array([+0.0000002462, -0.0002856080, +2.5831334692])}
        }
    }

bound_cuts = { 
    'H':{
        6: [455, 230],
        10: [250, 150],
        11: [600, 150],
        #13: [200, 600],
        13: [250, 750],
        14: [780, 150],
        16: [530, 200],
        17: [1000,100],
        20: [500, 150],
        21: [420, 290],
        22: [200, 150]},
    'K':{
        3:  [200, 1000],
        4:  [150, 250],
        13: [200, 400],
        14: [200, 400]}
}

class FitObjs:

    def __init__(self, s, x, u, continuum, watm_in, satm_in, mflux_in, 
                    mwave_in, mask, masterbeam, CRmask, initwave,
                    molmask):
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
        self.initwave = initwave
        self.molmask = molmask


class InParams:

    def __init__(self, inpath, outpath, initvsini, vsinivary, plotfigs, 
                    initguesses, bvcs, tagsA, tagsB, nights, mwave, mflux,
                    a0dict, xbounddict, maskdict):
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
        self.ips_tightmount_pars = ips_tight
        self.ips_loosemount_pars = ips_loose
        self.bound_cut_dic = bound_cuts

        
        self.methodvariance_tight = { 
        #                      6           14          16         21          22
            'H': np.array([0.0115497, 0.00834182, 0.00496473, 0.16538301, 0.00582415]),
        #                       4           5           6
            'K': np.array([0.00433336, 0.00113205, 0.00115613])
            }
        self.methodvariance_loose = { 
            'H': np.array([0.0115497, 0.00834182, 0.00496473, 0.16538301, 0.00582415]),
            'K': np.array([0.0225548, 0.00215941, 0.00119834])
            }



class InParamsA0:

    def __init__(self, inpath, outpath, plotfigs, tags, nights, humids,
                    temps, zds, press, obs, watm, satm, mwave, mflux,
                    cdbsloc, xbounddict, maskdict):
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
        self.ips_tightmount_pars = ips_tight
        self.ips_loosemount_pars = ips_loose
        self.bound_cut_dic = bound_cuts



class OrderDictCla:
    def __init__(self,):
        self.orderdict = { 
            'H':{
                1: [1.79350, 1.81560],
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
            'K':{
                1: [2.45317, 2.48280],
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

class TagStuffs:

    def __init__(self,night,watm,satm,a0contwave,continuum,ip):
        self.a0contwave = a0contwave
        self.continuum = continuum
        self.watm = watm
        self.satm = satm
        self.night = night
        self.ip = ip


def _setup_bound_cut(bound_cut_dic, band, order):
    """ Retrieve pixel bounds for where within each other significant telluric 
    absorption is present. If these bounds were not applied, analyzing some 
    orders would give garbage fits.

    Parameters
    ----------
    bound_cut_dic : Dict
        dict of pixel cuts on both sides (start & end) of the spectrum in 
        different orders
    band : str
        H or K band
    order : int
        spectrun order

    Returns
    -------
    list
        pixel cuts on both sides (start & end) of the spectrum in the given order
    """
    
    if band=='K':
        if int(order) in [3, 4, 13, 14]:
            bound_cut = bound_cut_dic[band][order]
        else:
            bound_cut = [150, 150]

    elif band=='H':
        if int(order) in [6, 10, 11, 13, 14, 16, 17, 20, 21, 22]:
            bound_cut = bound_cut_dic[band][order]
        else:
            bound_cut = [150, 150]
    
    return bound_cut
