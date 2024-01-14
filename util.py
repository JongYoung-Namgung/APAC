
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the APAC package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


import nibabel as nib
import numpy as np

def make_funcgii(dummy_file, input_arr, out_file):
	r"""Transfer array to func.gii file
	
	Parameters
	----------

	dummy_file: str 
		template func file, the dimension of the file must be same with input_arr
		func gifti file (extension: '.func.gii')

	input_arr: list | numpy.array
		input array to plot
		must be same dimension with dummy_file

	out_file: str
		file name for output func.gii file
		Extension ".func.gii" is no matter
	"""

	dummy = nib.load(dummy_file)
	dummy.darrays[0].data = input_arr.astype(np.float32)
	
	if '.'.join(out_file.split('.')[-2:]) == 'func.gii':
		pass
	else:
		out_file = out_file + '.func.gii'
	nib.save(dummy, out_file)
	print('Saved : {}'.format(out_file))