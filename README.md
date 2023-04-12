# Select ground motions using Weighted Multi-IM Objective


This code selects ground motions based on weighted multi-IM normalized objective function using different intensity measures including Sa, PGV, CAV, Ia, and D5-95


References: 
          
          Jawad Fayaz, Miguel Medalla, Pablo Torres, and Carmine Galasso (2023). "Recurrent Neural Networks based Generalized Ground Motion Model for Chilean Subduction Seismic Environment". Structural Safety, Vol. 100, Pages 102282


          Jawad Fayaz and Carmine Galasso (2022). "A Generalized Ground Motion Model for Consistent Mainshock-Aftershock Intensity Measures using Successive Recurrent Neural Networks". Bulletin of Earthquake Engineering, Vol. 20, Pages 6467-6486


INPUTS:

This code uses the target IMs provided in 'Target_IMs.xlsx' file to select best matching GMs using the IMs provided in 'IM_Data.mat' based on weighted multi-IM normalized objective function. The example files are provided with the code.


The 'IM_Data.mat' file which must be in the current folder and contain 2 variables

          'Spectra'       --> n x 1  Cell structure containing the spectra of 'n' GMs

          'Other_IMs'     --> (n+1) x 4  Cell structure containing the 4 IMs of 'n' GMs with 4 titles in the same order as the example file


The 'Target_IMs.xlsx' file must contain 3 columns in the same format as the example file. The column 1 can be empty, the code uses the columns 2 and 3 of the excel


OUTPUTS:

The outputs will be provided in 'SEL_IMs_GMs.mat' file which contain the scale and unscaled versions of the selected IMs. The output file will contain 

            'Sel_CAV_Unscaled','Sel_Spectra_Unscaled','Sel_Spectra_Scaled','Sel_PGV_Unscaled','Sel_Ia_Unscaled','Sel_CAV_Scaled','Sel_PGV_Scaled','Sel_Ia_Scaled','Sel_D595','Sel_Error','Sel_GMs_ids'


The variables are self-explanatory 
          
          'Sel_GMs_ids' --> n x 1  array containing indices of the selected GM IMs  with the same order as the 'IM_Data.mat' file
          
