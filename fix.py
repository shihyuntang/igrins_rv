from Engine.importmodule import *
import shutil

list_input = os.listdir('./Input')
list_input = [i for i in list_input if i[0] != '.']
list_input = [i for i in list_input if i[0] != '@']
list_input = [i for i in list_input if i != 'Prepdata']
list_input = [i for i in list_input if i != 'UseWv']

for tarname in list_input:
    nights = os.listdir(f'./Input/{tarname}')
    nights = [i for i in nights if i[0] != '.']
    nights = [i for i in nights if i[0] != '@']
    nights = [i for i in nights if i != 'std']
    for night in nights:
        list_A = os.listdir(f'./Input/{tarname}/{night}/A')
        for Afits in list_A:
            if 'sepc' in list_A:
                shutil.move(f'./Input/{tarname}/{night}/A/{Afits}',
                            f'./Input/{tarname}/{night}/A/{Afits[:18] + "spec" + Afits[-5:]}' )

        list_B = os.listdir(f'./Input/{tarname}/{night}/B')
        for Bfits in list_B:
            if 'sepc' in Bfits:
                shutil.move(f'./Input/{tarname}/{night}/B/{Bfits}',
                            f'./Input/{tarname}/{night}/B/{Bfits[:18] + "spec" + Bfits[-5:]}' )
