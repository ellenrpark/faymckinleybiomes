# faymckinleybiomes

GitHub repository to recreate the Fay and McKinley biomes from 2010-2024 (modified from Fay and McKinley, 2014; original biomes can be downloaded [here](https://doi.pangaea.de/10.1594/PANGAEA.828650)).

## Makig biomes
Run make_biomes.py to download and recreate these biomes.

**Notes:**
- User needs to specify output directory to download and store data files. If not, data will be downloaded in current working directory.
- Code creates directory /faymckinleybiomes/ where all data are downloaded. Below is the directory structure
```
├── faymckinleybiomes
|   make_biomes.py
|   ├── faymckinleybiomes
│     ├── data
│     ├── figures
│     ├── biomes_2010_<endyear>.nc
```
## Data sets
This workflow downloads and uses the following data from the following sources:
- Monthly 1ºx1º sea surface temperature (SST) and sea ice fraction from HadISST ([download here](https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html]), (Rayner et al., 2003)).
- Monthly 9km x 9km mapped MODIS-Aqua chlorophyll concentration from NASA's Ocean Color Level 3 & 4 browser ([download here](https://oceandata.sci.gsfc.nasa.gov/l3/)) downsampled to 1ºx1º
  - Note: To download these data an EarthData account and .netrc file are needed and must be properly configured on your device ([more info here](https://oceancolor.gsfc.nasa.gov/data/download_methods/))
- Argo 1ºx1º Mixed Layer climatology ([download here](https://mixedlayer.ucsd.edu/), (Holte et al., 2017))

## References
Fay, Amanda R; McKinley, Galen A (2014): Global Ocean Biomes: Mean and time-varying maps (NetCDF 7.8 MB) [dataset]. PANGAEA, https://doi.org/10.1594/PANGAEA.828650,
Supplement to: Fay, AR; McKinley, GA (2014): Global open-ocean biomes: mean and temporal variability. Earth System Science Data, 6(2), 273-284, https://doi.org/10.5194/essd-6-273-2014

Holte, J., L. D. Talley, J. Gilson, and D. Roemmich (2017), An Argo mixed layer climatology and database, Geophys. Res. Lett., 44, 5618–5626, doi:10.1002/2017GL073426.

Rayner, N. A.; Parker, D. E.; Horton, E. B.; Folland, C. K.; Alexander, L. V.; Rowell, D. P.; Kent, E. C.; Kaplan, A. (2003) Global analyses of sea surface temperature, sea ice, and night marine air temperature since the late nineteenth century J. Geophys. Res.Vol. 108, No. D14, 4407 10.1029/2002JD002670  (pdf ~9Mb)
