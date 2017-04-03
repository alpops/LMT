import scipy.io
from scipy.optimize import curve_fit
from numpy import histogram,ones,exp,savez
from fnmatch import fnmatch
from os import listdir
from optparse import OptionParser



def parse_obsnum_list(option, opt, value, parser):

	setattr(parser.values, option.dest, value.split(','))


def vlbi1mmFileSearch(obsnum, root='/data_lmt/vlbi1mm/', full=True):

	listDir = listdir(root)

	for ifile in listDir:
	
		if fnmatch(ifile, 'vlbi1mm*_%06d_*.nc' %obsnum):

			if full:

				ifile = root + ifile

			return ifile

	return ""

def form_full_variable(nc, data_name):

	data = nc.variables[data_name].data.copy()
	
	full_var = ''

	blank = ' '
	i = 0

	while (data[i] != blank):

		full_var += data[i]
		i+=1

	return full_var


def gauss(x,mu,sigma,A):

	return A*exp(-(x-mu)**2/2/sigma**2)


def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):

	return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)


def find_cal_ref(obsnum, nc_obs, max_scans=10, def_ref_A=5.0, def_ref_B=5.0):

	src_obs = form_full_variable(nc_obs, 'Header.Source.SourceName')
	ObsPgm = form_full_variable(nc_obs, 'Header.Dcs.ObsPgm')

	if ObsPgm != 'VlbiSched' and ObsPgm != 'Cal':

		print 'WARNING: %s scan %i is not of type VlbiSched or Cal' %(ObsPgm,obsnum) 

	ref_A = def_ref_A
	ref_B = def_ref_B

	is_cal = False

	j = obsnum
	calnum = j

	while ((not is_cal) and (abs(j-obsnum)<max_scans)):

		scan = vlbi1mmFileSearch(j)
		nc = scipy.io.netcdf_file(scan)
		
		ObsPgm_scan = form_full_variable(nc, 'Header.Dcs.ObsPgm')

		is_cal = (ObsPgm_scan == 'Cal')
	
		if(is_cal):

			calnum = j
			src = form_full_variable(nc, 'Header.Source.SourceName')

			print 'For %s scan %i (%s) found Cal scan %i (%s)' %(ObsPgm,obsnum,
										src_obs,calnum,src)

			if src != src_obs:
				print 'WARNING: Cal source is different from Observation source'

			ref_A,cold_A = fit_ref_voltage(nc.variables['Data.Vlbi1mmTpm.APower'].data.copy())
			ref_B,cold_B = fit_ref_voltage(nc.variables['Data.Vlbi1mmTpm.BPower'].data.copy())
			
		else:

			j -= 1

		nc.close()
				

	if (not is_cal):

		print 'ERROR: Cal scan not found within %i scans of Obs scan %i' %(max_scans,obsnum)
		print 'Reverting to default %i, %i ref voltages for A,B' %(def_ref_A, def_ref_B)
	
	return ref_A, cold_A, ref_B, cold_B, calnum


def fit_ref_voltage(cal_data, nbins=1000):

	vals,x  = histogram(cal_data, nbins)
	x = (x[1:]+x[:-1])/2.

	expected = (2., 1., 50., 6., 1., 50.)
	params, cov = curve_fit(bimodal, x, vals, expected)

	ref = max(params[0], params[3])
	cold = min(params[0], params[3])

	#print ref,cold
	return ref,cold


def calc_tsys(ref, sky_data, Tamb=280., cal=False):

	if cal:
		ref_data = ref
	else:
		ref_data = ref*ones(len(sky_data))
	
	Tsys = Tamb*(sky_data/(ref_data-sky_data))

	return Tsys



if __name__ == "__main__":


	parser = OptionParser()

	parser.add_option("-s", "--scans", type="string", action="callback", 
			callback=parse_obsnum_list, dest="scanlist",
			help="Input comma-separated list of observation scan numbers (obsnum)")
	
	#parser.add_option("-S", "--allscans", action="store", type="string", dest="scanbounds", 
	#		help="Input first and last observation scan number (obsnum)")

	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
			help="Prints time-averaged Tsys, source, and calibration obsnum\
				for each scan")

	parser.add_option("-o", "--output", action="store_true", dest="output",
			default=False, help="Outputs Tsys data into separate txt files")

	
	(options, args) = parser.parse_args()

	obsnum_list = map(int, options.scanlist)

	output = options.output
	verbose = options.verbose
			
	for obsnum in obsnum_list:

		obs_file = vlbi1mmFileSearch(obsnum)
        	nc_obs = scipy.io.netcdf_file(obs_file)

		src_obs = form_full_variable(nc_obs, 'Header.Source.SourceName')
 
		ref_A, cold_A, ref_B, cold_B, calnum = find_cal_ref(obsnum, nc_obs)

                if obsnum==calnum:
			
			cal_Tsys_A = calc_tsys(ref_A, cold_A, cal=True)
			cal_Tsys_B = calc_tsys(ref_B, cold_B, cal=True)
			
			if(verbose):

				print 'For Cal scan %i, Tsys_A = %f' %(calnum,cal_Tsys_A)
				print 'For Cal scan %i, Tsys_B = %f' %(calnum,cal_Tsys_B)
			continue 
			
			
		sky_A = nc_obs.variables['Data.Vlbi1mmTpm.APower'].data.copy()
		sky_B = nc_obs.variables['Data.Vlbi1mmTpm.BPower'].data.copy()

		t = nc_obs.variables['Data.Vlbi1mmTpm.Time'].data.copy()

		Tsys_A = calc_tsys(ref_A, sky_A)
		Tsys_B = calc_tsys(ref_B, sky_B)

		if(verbose):

			
			print 'Tsys A = %i K (time-avg)' %(Tsys_A.mean())
			print 'Tsys B = %i K (time-avg)' %(Tsys_B.mean()) 
			 

		if(output):

			savez('%i_tsys' %obsnum, Tsys_A = Tsys_A, Tsys_B = Tsys_B, t = t) 
				
		

