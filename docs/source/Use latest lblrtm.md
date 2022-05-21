# Use the Latest Version of lblrtm/lnfl/aer

*IGRINS RV Step 1 use `telfit` package to model atmosphere absolution lines (telluric lines), and `telfit` use the LBLTRM model to generate telluric spectrum. Currently, `telfit` is using lblrtm v12.2.1, lnfl v2.6, and aer_line_file v3.2 on [this zenodo](https://zenodo.org/record/1202479/#.YohJ7C-B2-Y). However, on May 20, 2022, the latest version are lblrtm v12.11, lnfl v3.2, and aer_line_file v3.8.1 on [Atmospheric and Environmental Research R&C (AER-RC)'s](https://github.com/AER-RC) Github repository.*

***

> Using the latest version of lblrtm/lnfl/aer **IS NOT REQUIRED** to run **IGRINS RV**, however, update lblrtm/lnfl/aer can help to fix several reported issue while running `telfit`.

## Benefits of using the latest version of lblrtm/lnfl/aer

1. Fix error while running Step 1 (calling `telfit` under the hood) using `ifort` compiler with lblrtm v12.2.1, lnfl v2.6, and aer_line_file v3.2.  
2. With the ability to run `telfit` with `ifort`, we can then fix the **End of File Runtime Error** (see issue [#9](https://github.com/kgullikson88/Telluric-Fitter/issues/9)). The **End of File Runtime Error** seems to occur when using `gfortran` compiler.


## Make sure the default installation of telfit is working

Before you try to update lblrtm/lnfl/aer files, please first make sure you followed {ref}`Packages installation (part 2) - Telfit` and that `telfit is working.

## Update telfit to use latest lblrtm/lnfl/aer

The lblrtm/lnfl/aer files `telfit` used are saved under a hidden folder `.Telfit` at the home dir.
Do 
```
~$ cd ~/.Telfit
```
to enter the folder. Under `./Telfit`, the folder structure is like:
```
~/.Telfit
├── aer_v_3.2.tar.gz
├── aerlbl_v12.2.1.tar.gz
├── aerlnfl_v2.6.tar.gz
├── lnfl
│   ├── build
│   ├── docs
│   └── ...
├── lblrtm
│   ├── build
│   ├── src
│   └── ...
├── aer_v_3.2
│   ├── line_files_By_Molecule
│   ├── line_file
│   └── ...
└── rundir#
└── ...
```

### 1. Getting the latest version of lblrtm/lnfl/aer

All files can be found on [Atmospheric and Environmental Research R&C (AER-RC)](https://github.com/AER-RC).\
More specifically, follow the clone/tarball_download instruction on each of the github README page below 
* lblrtm: [https://github.com/AER-RC/LBLRTM](https://github.com/AER-RC/LBLRTM)
* lnfl: [https://github.com/AER-RC/LNFL](https://github.com/AER-RC/LNFL)
* aer_line_file: [https://github.com/AER-RC/AER_Line_File](https://github.com/AER-RC/AER_Line_File)
  
to get folders: `LNFL-3.2`, `LBLRTM-12.11`, and `AER_Line_File` (take lblrtm v12.11, lnfl v3.2, and aer_line_file v3.8.1 for example).

> It is important to follow all instructions on their github README pages. Stuff related to `submodule` is critical to whether the update of lblrtm & lnfl is success or not!

### 2. Replacement with new version

For `lblrtm` and `lnfl`, simply delete all stuff under  `~/.Telfit/lblrtm/*` and `~/.Telfit/lnfl/*` then copy everything in the `LBLRTM-12.11/*` to `~/.Telfit/lblrtm/*` and everything in the `LNFL-3.2/*` to `~/.Telfit/lnfl/*`.

As for `AER_Line_File`, remember to follow the [Dependencies](https://github.com/AER-RC/AER_Line_File#dependencies) and [Downloading and Staging the Line File](https://github.com/AER-RC/AER_Line_File#downloading-and-staging-the-line-file) steps on [aer_line_file](https://github.com/AER-RC/AER_Line_File) to download the line files. 
If you follow the instruction, under `AER_Line_File` you will have:
```
./AER_Line_File/
├── LICENSE
├── get_line_file.py
├── aer_v_3.8.tar.gz
├── AER_Line_File
│   ├── RELEASE_NOTES_aer_linefile
│   ├── line_file
│   ├── spd_dep
│   ├── extra_brd_params
│   ├── line_files_By_Molecule
│   └── spectral_lines_for_MonoRTM
└── ...
```
Rename `./AER_Line_File/AER_Line_File` to `./AER_Line_File/aer_v_3.8`, then move the `aer_v_3.8` folder under `~/.Telfit/.`.

So, now, under your `~/.Telfit` you have:
```
~/.Telfit
├── aer_v_3.2.tar.gz
├── aerlbl_v12.2.1.tar.gz
├── aerlnfl_v2.6.tar.gz
├── lnfl
├── lblrtm
├── aer_v_3.2
├── aer_v_3.8
└── rundir#
└── ...
```

### 3. Update the setup.py

In `/Telluric-Fitter(-master)/setup.py` line [193](https://github.com/kgullikson88/Telluric-Fitter/blob/7ae98db278525e157d2d0abaf4697e2fe778d6bc/setup.py#L193), change
```python
linfile = u"{0:s}/aer_v_3.2/line_file/aer_v_3.2".format(TELLURICMODELING)
```
to
```python
linfile = u"{0:s}/aer_v_3.8/line_file/aer_v_3.8".format(TELLURICMODELING)
```
Then commented (#) out line [230](https://github.com/kgullikson88/Telluric-Fitter/blob/7ae98db278525e157d2d0abaf4697e2fe778d6bc/setup.py#L230) to line 237:
```python
    #Get/Unpack the tar files
    # for fname in ['aer_v_3.2.tar.gz', 'aerlnfl_v2.6.tar.gz', 'aerlbl_v12.2.1.tar.gz']:
    #     if fname not in os.listdir(TELLURICMODELING):
    #         url = '{}{}'.format(DATA_URL, fname)
    #         print('Downloading data from {} and putting it in directory {}'.format(url, TELLURICMODELING))
    #         download_file(url, '{}{}'.format(TELLURICMODELING, fname))
    #     print("Un-packing {}".format(fname))
    #     subprocess.check_call(["tar", "-xzf", '{}{}'.format(TELLURICMODELING, fname), '-C', TELLURICMODELING])
```
**so the next time we run `(igrins_rv) ~$ python setup.py build` the `./lnlf` and `./lblrtm` won't be overwritten by the old version!**

### 4. Re-run setup.py

Enter the `igrins_rv` environment (within which Telfit must be installed) and `cd` into `Telluric-Fitter(-master)`, then run
```
(igrins_rv) ~$ python setup.py build
(igrins_rv) ~$ python setup.py install
```

If you don't see any Error, then, congratulation, you successfully update your lblrtm/lnfl/aer versions!
