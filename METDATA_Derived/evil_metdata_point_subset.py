import netCDF4,numpy,sys,os
for fvar,ncvar in {'pr':'precipitation','rmean':'relative_humidity','srad':'surface_downwelling_shortwave_flux_in_air','tmmx':'air_temperature','tmmn':'air_temperature','vs':'wind_speed'}.items():
    mfd=netCDF4.MFDataset(os.path.join(sys.argv[1],fvar+'*.nc'),aggdim='day')
    for p in numpy.loadtxt(sys.argv[3]):
        numpy.savetxt(os.path.join(sys.argv[2],'%.3f_%.3f_%s.txt'%(p[1],p[0],fvar)),mfd.variables[ncvar][:,numpy.where(abs(p[1]-mfd.variables["lat"][:])<0.020833)[0][0],numpy.where(abs(p[0]-mfd.variables["lon"][:])<0.020833)[0][0]],fmt='%.3f',newline='\r\n',header='19790101',comments='')
