# frb_search

Package to create dynamic spectra from sky images (originally MWA) and then running de-dispersion with a specified DM range and creating de-dispersed time series (DTS).

Example usage:

 - create dynamic spectra from images in .fits format:
   create_dynaspec -w "(50,50)-(150,150)" -o ${gps} -f dirty_image_%dT%d_real.fits -v 10 -F -N 32 -X 0.100 -I 2 -C ${freq_ch}

   where gps is GPS time of the start of the MWA observation, -f provides the format of the input .fits file names, -C is the MWA coarse frequency channel.


 - de-disperse dynamic spectra:
 
   timeseries_all.sh "????_????" 0 0 5 "-S 1 -A -T -C 0.02893518518518518518 -I 2 -B ${freq_mhz} -E ${freq_ch} -P 0" 1 ${gps} 0 900 150

   where "????_????" is the naming of the directories with dynamic spectra, -C is the width of indiviual fine channel, -B provides center frequency -E provide channel number.

   The scripts calls program :

     dynaspec_search ${fits} ${out_fits}  -o ${gps} -l $dm_min -m $dm_max -s ${dm_step} -n ${snr_threshold} ${options}


If you use this software for your work, please cite the paper: https://ui.adsabs.harvard.edu/abs/2021PASA...38...26A/abstract
by including this citation into your bibliography file:

@ARTICLE{2021PASA...38...26A,
       author = {{Anderson}, G.~E. and {Hancock}, P.~J. and {Rowlinson}, A. and {Sokolowski}, M. and {Williams}, A. and {Tian}, J. and {Miller-Jones}, J.~C.~A. and {Hurley-Walker}, N. and {Bannister}, K.~W. and {Bell}, M.~E. and {James}, C.~W. and {Kaplan}, D.~L. and {Murphy}, Tara and {Tingay}, S.~J. and {Meyers}, B.~W. and {Johnston-Hollitt}, M. and {Wayth}, R.~B.},
        title = "{Murchison Widefield Array rapid-response observations of the short GRB 180805A}",
      journal = {\pasa},
     keywords = {gamma-ray bursts, gamma-ray bursts: individual: GRB 180805A, radio continuum: transients, neutron star mergers, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2021,
        month = jun,
       volume = {38},
          eid = {e026},
        pages = {e026},
          doi = {10.1017/pasa.2021.15},
archivePrefix = {arXiv},
       eprint = {2104.14758},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021PASA...38...26A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


