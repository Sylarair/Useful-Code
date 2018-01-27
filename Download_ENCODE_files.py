#!/usr/local/bin/python
# _*_ coding='utf-8' _*_
import os
import sys
import pandas as pd
import numpy as np
from pandas import DataFrame
import argparse
import math
import time
from multiprocessing import Pool
encode_data_path = '/data/projects/encode/data'
version = '0.1'
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputFile',help='the input information files for download datasets',required=True)
parser.add_argument('-o', '--outdir',help='the out directory for download datasets',required=True)
parser.add_argument('-a', '--assembly',help='Choosed genome',default='mm10')
parser.add_argument('-p', '--process',help='Number of processes',default='1')
args = parser.parse_args()

reference_gm = '/home/wangqin/reference/{}.chrom.sizes'.format(args.assembly)

def mkdir(path):
	if not os.path.exists(path):
		os.makedirs(path)
	else:
		pass

def Load_file(infile):
	os.system("grep -v 'bam' {} | grep -v 'unfiltered alignments' | grep -v 'signal p-value' | grep -v 'strand signal of all reads' > {}".format(infile, infile + '.temp.bed'))
	mat = pd.read_table(infile + '.temp.bed')#.values
	# mat = []
	return mat


# data = Load_file('/Users/wangqin/Downloads/e12.5day_control.tsv')
def rm_null_col(infile):
	data = Load_file(infile)

	for i in data.columns:
		if np.sum(data[i])!='nan':
			pass
		elif math.isnan(np.sum(data[i])):
			data.drop(data.drop(i),axis=1,inplace=True)

	return data

def download_files(infile):
	useful_col = ['File accession', 'File format', 'Output type', 'Experiment accession', 'Assay', 'Biosample term name', 'Biosample Age', 'Biosample treatments', 'Biosample subcellular fraction term name', 'Biosample phase', 'Biosample synchronization stage', 'Experiment target', 'Biological replicate(s)', 'File download URL', 'Assembly']

	print('{} starts to download data!'.format(time.ctime()))
	Array = Load_file(infile)
	Array1 = Array.values

	for i in range(len(Array1)):
		if Array['Assembly'][i] == args.assembly:
			id = Array['File accession'][i]
			postfix = Array['File format'][i]
			Type = Array['Output type'][i]
			exp_id = Array['Experiment accession'][i]
			exp_type = '_'.join(str(Array['Assay'][i]).split(' '))
			tissue = '_'.join(str(Array['Biosample term name'][i]).split(' '))
                        Date = '_'.join(str(Array['Biosample Age'][i]).split(' '))
                        treatment = '_'.join(str(Array['Biosample treatments'][i]).split(' '))
                        subcellular = '_'.join(str(Array['Biosample subcellular fraction term name'][i]).split(' '))
                        phase = '_'.join(str(Array['Biosample phase'][i]).split(' '))
                        synchronization_stage = '_'.join(str(Array['Biosample synchronization stage'][i]).split(' '))
                        target = Array['Experiment target'][i]
			rep = Array['Biological replicate(s)'][i]
			link = Array['File download URL'][i]

			keywords = [tissue, Date, treatment, subcellular, phase, synchronization_stage]

			if exp_type == 'ChIP-seq':
				exp_type = target.split('-')[0]
			else:
				exp_type = exp_type

			if postfix not in ['bam', 'bed', 'bigWig']:
				postfix = link.split('/')[-1].split('.')[-1]
                                if postfix == 'gz':
                                    postfix = 'bed.gz'
                                else:
                                    postfix = postfix
			else:
				postfix = postfix

			if rep == '1, 2' and Type == 'replicated peaks':
				rep = 'final'
			elif rep == '1, 2':
				rep = 'rep0'
			else:
				rep = 'rep' + rep

			tissue_last = '' + keywords[0]
                        for _itm in keywords[1:]:
				if str(_itm)=='nan':
					pass
				else:
					tissue_last += '_' + _itm

			store_data_path = args.outdir + '/' + exp_type + '/' + tissue_last
			mkdir(store_data_path)

			if os.path.exists('{}/{}/{}.{}'.format(encode_data_path, exp_id, id, postfix)):
				os.system('cp -r {} {}'.format('{}/{}/{}.{}'.format(encode_data_path, exp_id, id, postfix), store_data_path))
			else:
				os.system('wget -P {} -c {}'.format(store_data_path, link))

			os.system('cat {} > {}'.format('{}/{}.{}'.format(store_data_path, id, postfix), '{}/{}_{}.{}'.format(store_data_path, tissue_last, rep, postfix)))
		else:
			print('Line{} was skipped! It was {}'.format(i, Array['Assembly'][i]))

	print('{} ends to download data!'.format(time.ctime()))

def merge_bw(inpath, experiment, tis):
	prefix = '{0}/{1}/{2}/{2}'.format(inpath, experiment, tis)
	if os.path.exists('{}_rep0.bigWig'.format(prefix)):
		print('{}_rep0.bigWig already exists!'.format(prefix))
	else:
		rep_list_pre = os.popen("ls {}_rep[1-9].bigWig".format(prefix)).readlines()
		rep_list = [i.strip() for i in rep_list_pre]
		rep_num = len(rep_list)

		if rep_num == 0:
			print('{}/{}/{} does not exist bigWig files!'.format(inpath, experiment, tis))
		elif rep_num == 1:
			os.system('cat {} > {}'.format(rep_list[0], '{}_rep0.temp.bigWig'.format(prefix)))
		elif rep_num == 2:
			os.system('bigWigMerge {} {} {}'.format(rep_list[0], rep_list[1], '{}_rep0.biedGraph'.format(prefix)))
			os.system('sort -k1,1 -k2,2n {0}_rep0.bedGraph > {0}_rep0.sorted.bedGraph'.format(prefix))
			os.system('bedGraphToBigWig {0}_rep0.sorted.bedGraph {1} {0}_rep0.temp.bigWig'.format(prefix, reference_gm))
		elif rep_num == 3:
			os.system('bigWigMerge {} {} {} {}'.format(rep_list[0], rep_list[1], rep_list[2], '{}_rep0.biedGraph'.format(prefix)))
			os.system('sort -k1,1 -k2,2n {0}_rep0.bedGraph > {0}_rep0.sorted.bedGraph'.format(prefix))
			os.system('bedGraphToBigWig {0}_rep0.sorted.bedGraph {1} {0}_rep0.temp.bigWig'.format(prefix, reference_gm))
		else:
			print('{}/{}/{} exists more than 3 bigWig files!'.format(inpath, experiment, tis))

		os.system('rm -r {0}_rep0.bedGraph {0}_rep0.sorted.bedGraph'.format(prefix))

def merge_reps(args):
	exp_list = os.listdir(args.outdir)
	pool = Pool(processes=int(args.process))
	print('{} starts to merge data!'.format(time.ctime()))
	for exp in exp_list:
		tissue_list = os.listdir('{}/{}'.format(args.outdir, exp))
		for tissue in tissue_list:
			#merge_bw(args.outdir, exp, tissue)
			pool.apply_async(merge_bw, (args.outdir, exp, tissue, ))
			# prefix = '{0}/{1}/{2}/{2}'.format(args.outdir, exp, tissue)
			# if os.path.exists('{}_rep0.bigWig'.format(prefix)):
			# 	print('{}_rep0.bigWig already exists!'.format(prefix))
			# else:
			# 	rep_list_pre = os.popen("ls {}_rep[1-9].bigWig".format(prefix)).readlines()
			# 	rep_list = [i.strip() for i in rep_list_pre]
			# 	rep_num = len(rep_list)
            #
			# 	if rep_num == 0:
			# 		print('{}/{}/{} does not exist bigWig files!'.format(args.outdir, exp, tissue))
			# 	elif rep_num == 1:
			# 		os.system('cat {} > {}'.format(rep_list[0], '{}_rep0.temp.bigWig'.format(prefix)))
			# 	elif rep_num == 2:
			# 		os.system('bigWigMerge {} {} {}'.format(rep_list[0], rep_list[1], '{}_rep0.biedGraph'.format(prefix)))
			# 		os.system('sort -k1,1 -k2,2n {0}_rep0.bedGraph > {0}_rep0.sorted.bedGraph'.format(prefix))
			# 		os.system('bedGraphToBigWig {0}_rep0.sorted.bedGraph {1} {0}_rep0.temp.bigWig'.format(prefix, reference_gm))
			# 	elif rep_num == 3:
			# 		os.system('bigWigMerge {} {} {} {}'.format(rep_list[0], rep_list[1], rep_list[2], '{}_rep0.biedGraph'.format(prefix)))
			# 		os.system('sort -k1,1 -k2,2n {0}_rep0.bedGraph > {0}_rep0.sorted.bedGraph'.format(prefix))
			# 		os.system('bedGraphToBigWig {0}_rep0.sorted.bedGraph {1} {0}_rep0.temp.bigWig'.format(prefix, reference_gm))
			# 	else:
			# 		print('{}/{}/{} exists more than 3 bigWig files!'.format(args.outdir, exp, tissue))
            #
			# 	os.system('rm -r {0}_rep0.bedGraph {0}_rep0.sorted.bedGraph'.format(prefix))
	pool.close()
	pool.join()

	print('{} ends to merge data!'.format(time.ctime()))

def main():
	info_file = args.inputFile #'/Users/wangqin/Downloads/e12.5day_control.tsv'
	# download_files(info_file)
	merge_reps(args)

if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		print('Sorry! Keyboard interrupts me! Exit!')
		sys.exit(0)
