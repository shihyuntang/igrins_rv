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
    for night in nights:
        list_A = os.listdir(f'./Input/{tarname}/{night}/A')
        for Afits in list_A:
            if 'sepc' in list_A:
                shutil.move(Afits, Afits[:18] + 'spec' + Afits[-5:])

        list_B = os.listdir(f'./Input/{tarname}/{night}/B')
        for Bfits in list_B:
            if 'sepc' in Bfits:
                shutil.move(Bfits, Bfits[:18] + 'spec' + Bfits[-5:])
