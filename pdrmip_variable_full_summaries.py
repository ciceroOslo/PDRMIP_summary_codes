## Script to produce variable summaries for pdrmip
## This produces full summaries (means and standard deviations)
## for each of a selected list of variables
## for each model and experiment
## Three csv files are produced for each variable
## containing tabulated csv-data of the mean 
## the standarad deviation and a formatted version with both 
## mean and standard deviation tabulated by experiment and model.
##
## Under version 2, the cdo module has been removed. 
## To run this script only the python module needs to be loaded.



## Importing packages
import sys, os, glob, netCDF4, csv, math, pickle
import numpy as np
from operator import itemgetter

#Selected variables in different categories:
vars_Energy_fluxes = ['hfls_Amon', 'hfss_Amon', 'rlds_Amon', 'rlus_Amon', 'rsus_Amon', 'rsds_Amon', 'rsdt_Amon', 'rsut_Amon', 'rlut_Amon']

vars_Clear_sky_rad_Fluxes = ['rsutcs_Amon', 'rsdscs_Amon', 'rsuscs_Amon', 'rlutcs_Amon', 'rldscs_Amon']#'rluscs_Amon']

vars_meterological_fields = ['tas_Amon', 'ts_Amon', 'pr_Amon', 'prc_Amon', 'uas_Amon', 'vas_Amon', 'clt_Amon', 'ps_Amon', 'huss_Amon', 'hurs_Amon', 'evspsbl_Amon', 'prw_Amon']

vars_extra = ['abs550aer_Amon', "od550aer_Amon", "loadbc", "loadso4"]

vars = vars_Energy_fluxes + vars_Clear_sky_rad_Fluxes + vars_meterological_fields + vars_extra

#Ocean varieties:
ocean = ['coupled', 'fsst']

#All experiments:
experiments = ['base', 'base2','co2x2', 'ch4x3', 'solar', 'bcx10', 'sulx5', 'bcx10asia', 'sulx10asia', 'sulx10eur', 'sulred', 'sulasiared', 'cfc12', 'cfc11', 'n2o1p', 'ozone', 'lndus', 'bcslt'] 
# ,'base3', 'bcsl2']

#Just core experiments
#experiments = ['base', 'co2x2', 'ch4x3', 'solar', 'bcx10', 'sulx5']
#experiments = ['base', 'co2x2']
#ozone runs:
#experiments = ['ozone', 'o3all']

#Phase 2 experiments
#experiments = ['base2', 'cfc12', 'cfc12', 'cfc11', 'n2o1p', 'ozone', 'lndus', 'bcslt']

#Empty list to put the models
models = list()

# Reading model names from the folder names in renamed:
# This is where the files are at our local system
# To run somewhere else you will have to adjust the paths
# Also the organisation from this is folders per model
# then folders per ocean varieties
# then all files under there
# If you want to run this code elsewhere, we suggest you
# organise the data in the same way, or set up a symlink tree
# that is organised this way to do the calculations as was
# done in the nird repository
root, dirs, files = next(os.walk("/div/pdo/pdrmip/renamed/renamed/"))#.next()

# Adding models to list by folder names
for name in dirs:
    models.append(name)

## Select out specific models (for testing)
#models = itemgetter(1,9)(models)

models = sorted(models)

print(models)

#Some models have differing conventions:
#Make array of factors to account for this:
"""
fix_factor = np.ones((len(vars),len(models)))
#clt as fractions rather than percentages:
fix_factor[vars.index('clt_Amon'), models.index('CanESM2')] = 100
fix_factor[vars.index('clt_Amon'), models.index('HadGEM2')] = 100
fix_factor[vars.index('clt_Amon'), models.index('HadGEM3')] = 100
#evspsbl opposite sign:
fix_factor[vars.index('evspsbl_Amon'), models.index('CanESM2')] = -1
# hfss opposite sign
fix_factor[vars.index('hfss_Amon'), models.index('ECHAM-HAM')] = -1
fix_factor[vars.index('hfss_Amon'), models.index('HadGEM3')] = -1
"""
## Creating a mean data vector in which to store the calculated values:
meandata = np.zeros((len(vars),len(models), len(experiments), len(ocean)))
diffdata = np.zeros((len(vars),len(models), len(experiments), len(ocean)))
sigmadata = np.zeros((len(vars),len(models), len(experiments), len(ocean)))
print("number of variables, models, experiments, oceans: " + str(meandata.shape)) 

### Making headers for tables to be printed:
## Experiment header needs to be a formatted string
## with each experiment comma separated:
header_exp = "Models"
for e in experiments:
    header_exp = header_exp + "," + e

print(header_exp)
print(models)

## The model header needs to be a string coloumn vector:
models_header = np.array(models)[:,np.newaxis]
print(models_header)

## Process only mean values
mean_only = False

# Defining the function area weighted mean calculation
def area_mean (Array, mod_lon, mod_lat, yms, yme):
    print(Array.shape)
    Array[Array > 1e10] = 0 # Change any values exponentially high to zero, to avoid missing values messing up the means.
    if(len(Array.shape) == 4): # Check if the array is 4-dimentional.
        Array3D = Array.mean(1) # mean over the altitudes for 4D arrays. Sometimes there is only one layer, but the netcdf exist as 4 dimentions.
    else:
        Array3D = Array # For 3-dimention arrays
    tdim = Array3D.shape[0] #Gives the length of the first dimention (time)
    tseq = range(tdim)[yms:yme+1]
    tgl = [0] * len(tseq)
    n = 0
    for i in tseq:
        mat = Array3D[i,:,:]
        nlon=mat.shape[1]
        nlat=mat.shape[0]
        middellat = [0] * nlat
        for ilat in range(nlat):
            for ilon in range(nlon):
                middellat[ilat] = middellat[ilat] + mat[ilat, ilon]/nlon
        sumcos = 0
        gl = 0
        for ilat in range(nlat):
            fi = mod_lat[ilat]*math.pi/180
            gl = gl + middellat[ilat]*math.cos(fi)
            sumcos = sumcos + math.cos(fi)
        tgl[n] = tgl[n] + gl/sumcos
        n = n + 1
    return tgl;


def area_mean_fast (Array, mod_lon, mod_lat, yms, yme):
    print("There are " + str(len(Array.shape)) +  " dimentions with following lengths:")
    print(Array.shape)
    Array[Array > 1e10] = 0 # Change any values exponentially high to zero, to avoid missing values messing up the means. 
    if(len(Array.shape) == 4): # Check if the array is 4-dimentional.
        Array3D = Array.mean(1) # mean over the altitudes for 4D arrays. Sometimes there is only one layer, but the netcdf exist as 4 dimentions.
    else:
        Array3D = Array # For 3-dimention arrays
    
    tdim = Array3D.shape[0] #Gives the length of the first dimention (time)
    tseq = range(tdim)[yms:yme+1] #Picks out the specified indices in time dimention specified by the start and end year
  
    Array3D2 = Array3D[tseq,:,:] #Picks out the specified matries from the selected indices into new array
    mat = Array3D2.mean(0) #Time mean over the new array
    tgl = [0] * len(tseq) #Creates a vector of zeros with the length of the tseq 
    
    nlon=mat.shape[1] # Gives the number of longitude values
    nlat=mat.shape[0] # Gives the number of latitude values
    meanlat = [0] * nlat # Creates a vector of zeros with the length = nlat
    for ilat in range(nlat):   # iteration over latitude sequence
        for ilon in range(nlon): # iteration over longitude sequence
            meanlat[ilat] = meanlat[ilat] + mat[ilat, ilon]/nlon # averaging the latitude over the longitudes  
    sumcos = 0
    gl = 0
    for ilat in range(nlat): 
        fi = math.cos(mod_lat[ilat]*math.pi/180) # Calculating the weighting factor
        gl = gl + meanlat[ilat]*fi # Field weighting mean latitudes  
        sumcos = sumcos + fi

    tgl = gl/sumcos # Field mean
    return tgl;

list_several_files = []
list_missing_data = []

## Looping over variables, ocean, models and experiments
for v in vars:
    for o in ocean:
        for m in models:
            print("")
            print("MODEL: "+ m)
            for e in experiments:
                #Getting the folder name and fetching the matching files
                folder = '/div/pdo/pdrmip/renamed/renamed/%s/%s/'%(m,o)
                files = glob.glob( '/div/pdo/pdrmip/renamed/renamed/%s/%s/%s_*pdrmip-%s_*.nc'%(m,o,v,e))

                #Skipping faulty/outdated files
                if m == 'ECHAM-HAM' and o == 'coupled' and v == 'rlus' and e == 'bcsl2':
                    continue

                # Handling situations with extra files
                # This part may need to be changed a bit
                # if symlinks are used
                if(len(files) > 1):
                    links = 0
                    for f in files:
                        if os.path.islink(f):
                            links = links +1
                            files.remove(f)
                        elif f.endswith('T42.nc'):
                            print('T42 regridded file %s'%f)
                            files.remove(f)  
                    if links == 0:
                        list_several_files.append("%s_%s_%s_%s"%(v,e,m,o))
                        if o == "coupled" and m == "NorESM1" and e == "bcslt":
                            print("NorESM1 acting up as expected")
                            f = files[0]
                            f2 = files[1]
                        elif len(files)<2:
                            print("Superfluous file removed")
                            f = files[0]
                        else:
                             print("There are too many files")
                             print(files)
                             print("%s %s %s %s "%(v,e,m,o))
                             sys.exit(4)
                    elif links == 1:
                        f = files[0]
                        print("No problem, just symlink")
                    else:
                        print("Something wrong! Only symlinks for %s %s %s %s"%(v,e,m,o))              
                        sys.exit(4)
                elif(len(files) < 1):
                    list_missing_data.append("%s_%s_pdrmip-%s_%s"%(v,m,e,o))
                    continue
                    
                else:
                    f = files[0]
                # Finding the year span for which the data is defined:
                # The try except takes care of some files on the wrong 
                # format that shouldn't really be in the folder, but are
                try: 
                    yearsStart = int(f.split("_")[-1][0:4])
                    yearsEnd = int(f.split("_")[-1][7:11])
                except ValueError:
                    print('Something wrong with file ' + f)
                    continue 
                        
                # Depending on whether the run is coupled or fsst,
                # we choose a different subset of years to mean over:
                if(o == 'coupled'):
                    ## For Coupled:
                    ## add 50 years to start year  
                    ## add 99 years and 11 months to get end year
                    ysadd = 50
                    yeadd = 99
                    ys = yearsStart + ysadd
                    ye = yearsStart + yeadd
                    yms = ysadd*12
                    yme = (yeadd*12)+11
                    ## Special case for NorESM1 experiment bcslt:
                    ## Due to incomplete model run, the mean is calculated 
                    ## for 32 months earlier than other models.
                    if(m == "NORESM1" and e == "bcslt"):
                        yme = 66*12 + 11
                        yms2 = 0
                        yme2 = 20*12 + 11
                else:
                    ## For fSST:
                    ## add 5 years to start year  
                    ## add 14 years and 11 months to end year  
                    ysadd = 5
                    yeadd = 14
                    ys = yearsStart + ysadd
                    ye = yearsStart + yeadd
                    yms = ysadd*12
                    yme = (yeadd*12)+11

                print(ys)
                print(ye)
                print(yms)
                print(yme)
                print(f)
                        
                # At this point the extraction from the file is 
                # processed by the netcdf package
                # read in the netcdf file
                nc = netCDF4.Dataset(f,'r')

                #Split off the Amon/day part of the variable name
                var_name = v.split("_")[0]
                print(var_name)
                # Retrive variable values from file
                nc_var = nc.variables[var_name]
                #print(nc_var)
                var_array = np.array(nc_var)
                        
                # Retrive lon and lat values coordinates from file
                nc_lon = nc.variables['lon']
                nc_lat = nc.variables['lat']
                lon_array = np.array(nc_lon)
                lat_array = np.array(nc_lat)                    
                        
                nc.close()
                if(m == "NORESM1" and e == "bcslt" and o == "coupled"):
                    # read in the second netcdf file
                    nc = netCDF4.Dataset(f2,'r')
                
                    # Retrive variable values from file
                    nc_var = nc.variables[var_name]
                    #print(nc_var)
                    var_array2 = np.array(nc_var)
                        
                    # Retrive lon and lat values coordinates from file
                    nc_lon = nc.variables['lon']
                    nc_lat = nc.variables['lat']
                    lon_array2 = np.array(nc_lon)
                    lat_array2 = np.array(nc_lat)                    
                        
                    nc.close()
                    
                if(mean_only):
                    # Spatial and temporal mean
                    mean = area_mean_fast(var_array, lon_array, lat_array, yms, yme)
                    sigma = np.nan
                    #Add results from other file for NorESM1, coupled bcslt
                    if(m == "NORESM1" and e == "bcslt" and o == "coupled"):
                        mean = mean + area_mean_fast(var_array2, lon_array2, lat_array2, yms2, yme2) 
                else:
                    # Spatial mean:
                    fldmean = area_mean(var_array, lon_array, lat_array, yms, yme)
                    #Add results from other file for NorESM1, coupled bcslt
                    if(m == "NORESM1" and e == "bcslt" and o == "coupled"):
                        fldmean = fldmean + area_mean(var_array2, lon_array2, lat_array2, yms2, yme2) 
                    # Save field mean time series to pickle file
                    pkl_filename = 'fldmeans/' + v + '_' +  m + '_pdrmip-' + e + '_' + o + '_' + str(ys) + '-' + str(ye) + '.pkl'
                    print(pkl_filename)
                    pkl_fileObject = open(pkl_filename, 'wb')
                    pickle.dump(fldmean, pkl_fileObject)
                    pkl_fileObject.close()
                    # Time mean:
                    mean = np.mean(fldmean)
                    # Time standard deviation
                    sigma = np.std(fldmean)
          
                print("Mean: " + str(mean))
                print("Stdev: "+ str(sigma))
                # Then we store the mean value in the appropriate place in the matrix:
                meandata[vars.index(v), models.index(m), experiments.index(e), ocean.index(o)] = mean
                sigmadata[vars.index(v), models.index(m), experiments.index(e), ocean.index(o)] = sigma

                ## Then store the data with only differences from base
                #if(e is 'base'):
                #    diffdata[vars.index(v), models.index(m), experiments.index(e), ocean.index(o)] = mean
                #else:
                #    diffdata[vars.index(v), models.index(m), experiments.index(e), ocean.index(o)] = mean -  meandata[vars.index(v), models.index(m), experiments.index('base'), ocean.index(o)]
                                                
                #Printing filename for symlinks
            """
            else:
                print(f)
  
            """
       # Just printing tables:
      
        #Subsetting a 2-dimensional part of the matrix containing mean data for a particular variable
        #and ocean setup
        towrite = meandata[vars.index(v), :,:, ocean.index(o)].reshape(len(models), len(experiments))
        print(header_exp)
        print(towrite)
        print(models_header)
        print(np.hstack((towrite, models_header)))
        #Making an appropriate filename to write the table to:
        filename = './results/masan_results/%s_%s_meandata.csv'%(v, o)

        # Then saving the mean data for this variable and ocean type with headers
        np.savetxt(filename, np.hstack((models_header, towrite)), fmt = '%s', delimiter = ',', header = header_exp , comments = '')

        #Subsetting a 2-dimensional part of the matrix containing diff data for a particular variable
        #and ocean setup
        #towrite = diffdata[vars.index(v), :,:, ocean.index(o)].reshape(len(models), len(experiments))
        
        #Making an appropriate filename to write the table to:
        #filename = './results/masan_results/%s_%s_diffdata.csv'%(v, o)

        # Then saving the difference data for this variable and ocean type with headers
        #np.savetxt(filename, np.hstack((models_header, towrite)), fmt = '%s', delimiter = ',', header = header_exp)
   
        #Subsetting a 2-dimensional part of the matrix containing standard deviation data
        #for a particular variable and ocean setup
        towrite = sigmadata[vars.index(v), :,:, ocean.index(o)].reshape(len(models), len(experiments))
        
        #Making an appropriate filename to write the table to:
        # Particular place we put results, change this for
        # your system
        filename = './results/masan_results/%s_%s_sigmadata.csv'%(v, o)

        # Then saving the standard deviation data for this variable and ocean type with headers
        np.savetxt(filename, np.hstack((models_header, towrite)), fmt = '%s', delimiter = ',', header = header_exp, comments = '')     

        # Particular place we put results, change this for
        # your system
        filename_nice = "./results/masan_results/%s_%s_meanwstd.csv"%(v,o)
        printout_nice = [[header_exp]]
        print(printout_nice)
        for m in models:
            row = [m]
            for e in experiments:
                row.append("%.6G \u00B1 %.3G"%(meandata[vars.index(v), models.index(m), experiments.index(e), ocean.index(o)], sigmadata[vars.index(v), models.index(m), experiments.index(e), ocean.index(o)]))
            printout_nice.append(row)
        with open(filename_nice, 'w') as fi:
            writer = csv.writer(fi, delimiter = ',')
            writer.writerows(printout_nice)
print(list_several_files)

# Also writing out some files that keep track of
# missing data, and data with several files for
# debugging
with open("missing_and_metadata.txt", 'w') as fi:
    fi.write("Missing data for: ")
    for me in list_missing_data:
        fi.write("\n %s"%me)
    fi.write("More than one file: ")
    for ef in list_several_files:
        fi.write("\n %s"%ef)
