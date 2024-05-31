R scripts used for data analysis in the manuscript: 'Extreme compound climate hazards surge by 2040: governance is critical for the world’s most vulnerable'

by Craparo, A.C.W., Minoarivelo, O. H., Basel, A.M., Nguyen, K.T., Dao., H., Birner, J., Yonetani, M., Antinoja, E.S., Sanchez Torres, D.G., Brown, O.D., Pacillo., G., Harper, A., Läderach., P.

--------------------------------------------------------
The codes were run through RStudio 2023.06.1 using R version 4.3.1 on a device with a Windows 11 64-bit OS, an Intel Core i7-10510U @1.80GHz 2.30 GHz processor and a 16Gb RAM.
The processing time for all codes to run is about 30 minutes.

All raw input data (climate, FDSP population data and location) are available in the 'input_data' folder.

In order to replicate the analytics of the manuscript, the codes should be executed in the following order:

(1) The script named 'individual_hazard_treatment.R' creates the individual hazards (heat, drought, flood) indices. Figures output from the script are: Fig. 4 (B, C, D), Fig. 5 (B, C, D), Fig. S1, Fig. S2, and Fig. S4.

(2) The script named 'composite_hazard_build.R' creates the composite index. It needs to be executed after 'individual_hazard_treatment.R'. Figure output from the script is Fig. S3.

(3) The script named 'FDSP_exposure_compoud.R' is doing the statistical analyses with regard to the number of FDSP within each class of compound climate hazard, both with and without governance. It needs to be excecuted after 'composite_hazard_build.R'. Figures output from the script are: Fig. 1, Fig. 2, Fig. 3.

(4) The scrip named 'FDSP_exposre_individual_hazards.R' is doing the statistical analyses with regard to the number of FDSP within each class of individual climate hazards (heat, drought and flood). It needs to be excecuted after 'composite_hazard_build.R'. Figures output from this script are: Fig. 4A and Fig. 5A.
(5) To generate outputs corresponding to SSP126, the codes with the suffix 'ssp126' should be run following the same order as above. These codes are located in the 'ssp126' branch in here.
