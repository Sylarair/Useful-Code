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

encode_data_path = '/data/projects/encode/data'
version = '0.1'
parser = argparse.ArguementParser()
parser.add_arguement('-i', '--inputFile',description='the input information files for download datasets',required=True)
parser.add_arguement('-O', '--outdir',description='the out directory for download datasets',required=True)
parser.add_arguement('-A', '--assembly',description='Choosed genome',default='mm10')
args = parser.parse_args()

def mkdir(path):
	if not os.path.exists(path):
		os.makedirs(path)
	else:
		pass

def Load_file(infile):
	mat = pd.read_table(infile, nrows=100)#.values
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
			exp_type = Array['Assay'][i]
			tissue = Array['Biosample term name'][i]
			Date = Array['Biosample Age'][i]
			treatment = Array['Biosample treatments'][i]
			subcellular = Array['Biosample subcellular fraction term name'][i]
			phase = Array['Biosample phase'][i]
			synchronization_stage = Array['Biosample synchronization stage'][i]
			target = Array['Experiment target'][i]
			rep = Array['Biological replicate(s)'][i]
			link = Array['File download URL'][i]

			keywords = [tissue, Date, treatment, subcellular, phase, synchronization_stage]

			if exp_type == 'ChIP-seq':
				exp_type = target.split('-')[0]
			else
				exp_type = exp_type

			if postfix not in ['bam', 'bed', 'bigWig']:
				postfix = link.split('/').split('.')[-1]
			else:
				postfix = postfix

			if rep == '1, 2' and Type == 'replicated peaks':
				rep = 'final'
			elif rep == '1, 2'
				rep = 'rep0'
			else:
				rep = 'rep' + rep

			tissue_last = ''
			for _itm in keywords:
				if math.isnan(_itm):
					pass
				else:
					tissue_last += '_' + itm

			store_data_path = args.outdir + '/' + exp_type + '/' + tissue_last
			mkdir(store_data_path)

			if os.path.exists('{}/{}/{}.{}'.format(encode_data_path, exp_id, id)):
				os.system('cp -r {} {}'.format('{}/{}/{}.{}'.format(encode_data_path, exp_id, id, postfix), store_data_path))
			else:
				os.system('wget -P {} -c {}'.format(store_data_path, link))

			os.system('cat {} > {}'.format('{}/{}.{}'.format(store_data_path, exp_id, postfix), '{}/{}_{}.{}'.format(store_data_path, tissue_last, rep, postfix)))
		else:
			print('Line{} was skipped! It was {}'.format(i, Array['Assembly'][i]))

	print('{} ends to download data!'.format(time.ctime()))

def main():
	info_file = args.inputFile #'/Users/wangqin/Downloads/e12.5day_control.tsv'
	download_files(info_file)

if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		print('Sorry! Keyboard interrupts me! Exit!')
		sys.exit(0)












