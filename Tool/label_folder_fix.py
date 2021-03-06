from Engine.importmodule import *
import shutil

print('Program Start..')


list_input = os.listdir('./Input')
list_input = [i for i in list_input if i[0] != '.']
list_input = [i for i in list_input if i[0] != '@']
list_input = [i for i in list_input if i != 'Prepdata']
list_input = [i for i in list_input if i != 'UseWv']
list_input = [i for i in list_input if i != 'TauBoo']

print(f'Process {list_input} folders')
time.sleep(5)

for target_folder in list_input:
    print(f'Doing {target_folder}')

    master_log    = pd.read_csv('./Engine/IGRINS_MASTERLOG.csv')
    star_files    = master_log[(master_log['OBJNAME'].str.contains(target_folder, regex=True, na=False)) &
                               (master_log['OBJTYPE'].str.contains('TAR',         regex=True, na=False)) ]

    star_files    = Table.from_pandas(star_files)


    n = 1
    while len(star_files['CIVIL']) == 0:
        starnew = target_folder[:n]+' '+target_folder[n:]
        star_files = master_log[(master_log['OBJNAME'].str.contains(starnew, regex=True, na=False)) &
                                (master_log['OBJTYPE'].str.contains('TAR',   regex=True, na=False)) ]
        n += 1
        print(starnew)
        if n == len(target_folder)+1:
            sys.exit('TARGET NAME NOT FOUND IN CATALOG - CHECK INPUT!')

    print(f'We have {star_files["CIVIL"]}')

    for dateUT in np.unique(star_files['CIVIL']):
        print(f'Doing {dateUT}')
        date_star_files = star_files[ star_files['CIVIL'] == dateUT ]

        for beam, tag in zip(date_star_files['FRAMETYPE'], date_star_files['FILENUMBER']) :
            # Check A folder...
            if beam == 'A':
                if not os.path.isdir(f'./Input/{target_folder}/{dateUT}/A/SDCH_{dateUT}_{int(tag):04d}.sn.fits'):
                    try:
                        shutil.move(f'./Input/{target_folder}/{dateUT}/B/SDCH_{dateUT}_{int(tag):04d}.sn.fits',
                                    f'./Input/{target_folder}/{dateUT}/A/SDCH_{dateUT}_{int(tag):04d}.sn.fits')
                        shutil.move(f'./Input/{target_folder}/{dateUT}/B/SDCH_{dateUT}_{int(tag):04d}.spec.fits',
                                    f'./Input/{target_folder}/{dateUT}/A/SDCH_{dateUT}_{int(tag):04d}.spec.fits')
                        shutil.move(f'./Input/{target_folder}/{dateUT}/B/SDCK_{dateUT}_{int(tag):04d}.sn.fits',
                                    f'./Input/{target_folder}/{dateUT}/A/SDCK_{dateUT}_{int(tag):04d}.sn.fits')
                        shutil.move(f'./Input/{target_folder}/{dateUT}/B/SDCK_{dateUT}_{int(tag):04d}.spec.fits',
                                    f'./Input/{target_folder}/{dateUT}/A/SDCK_{dateUT}_{int(tag):04d}.spec.fits')
                        print(f'Move B/SDCK_{dateUT}_{int(tag):04d} to A/SDCK_{dateUT}_{int(tag):04d}')
                    except:
                        print(f'     --> NO FILE ./Input/{target_folder}/{dateUT}/B/SDCH_{dateUT}_{int(tag):04d}')
                else:
                    pass
            elif beam == 'B':
                if not os.path.isdir(f'./Input/{target_folder}/{dateUT}/B/SDCH_{dateUT}_{int(tag):04d}.sn.fits'):
                    try:
                        shutil.move(f'./Input/{target_folder}/{dateUT}/A/SDCH_{dateUT}_{int(tag):04d}.sn.fits',
                                    f'./Input/{target_folder}/{dateUT}/B/SDCH_{dateUT}_{int(tag):04d}.sn.fits')
                        shutil.move(f'./Input/{target_folder}/{dateUT}/A/SDCH_{dateUT}_{int(tag):04d}.spec.fits',
                                    f'./Input/{target_folder}/{dateUT}/B/SDCH_{dateUT}_{int(tag):04d}.spec.fits')
                        shutil.move(f'./Input/{target_folder}/{dateUT}/A/SDCK_{dateUT}_{int(tag):04d}.sn.fits',
                                    f'./Input/{target_folder}/{dateUT}/B/SDCK_{dateUT}_{int(tag):04d}.sn.fits')
                        shutil.move(f'./Input/{target_folder}/{dateUT}/A/SDCK_{dateUT}_{int(tag):04d}.spec.fits',
                                    f'./Input/{target_folder}/{dateUT}/B/SDCK_{dateUT}_{int(tag):04d}.spec.fits')
                        print(f'Move A/SDCK_{dateUT}_{int(tag):04d} to B/SDCK_{dateUT}_{int(tag):04d}')
                    except:
                        print(f'     --> NO FILE ./Input/{target_folder}/{dateUT}/A/SDCH_{dateUT}_{int(tag):04d}')
                else:
                    pass

print('Done')
