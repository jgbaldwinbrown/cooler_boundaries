#!/usr/bin/env python3

# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import python package for working with cooler files and tools for analysis
import cooler
import cooltools.lib.plotting
from cooltools import insulation

from packaging import version
if version.parse(cooltools.__version__) < version.parse('0.5.4'):
	raise AssertionError("tutorials rely on cooltools version 0.5.4 or higher,"+
						 "please check your cooltools version and update to the latest")

import cooltools
import itertools
from matplotlib.ticker import EngFormatter

from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bioframe
from skimage.filters import threshold_li, threshold_otsu
import bbi

def getdata():
	data_dir = './data/'
	cool_file = cooltools.download_data("HFF_MicroC", cache=True, data_dir=data_dir)

	resolution = 10000
	clr = cooler.Cooler(f'{data_dir}test.mcool::resolutions/{resolution}')
	windows = [3*resolution, 5*resolution, 10*resolution, 25*resolution]
	insulation_table = insulation(clr, windows, verbose=True)

	norm = LogNorm(vmax=0.1, vmin=0.001)
	return data_dir, resolution, clr, windows, insulation_table, norm

def print_first_summary(windows, insulation_table):
	first_window_summary =insulation_table.columns[[ str(windows[-1]) in i for i in insulation_table.columns]]
	print(insulation_table[['chrom','start','end','region','is_bad_bin']+list(first_window_summary)].iloc[1000:1005])

# Functions to help with plotting
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
	start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
	n = matrix_c.shape[0]
	t = np.array([[1, 0.5], [-1, 0.5]])
	matrix_a = np.dot(np.array([
		(i[1], i[0]) for i in itertools.product(start_pos_vector[::-1], start_pos_vector)
	]), t)
	x = matrix_a[:, 1].reshape(n + 1, n + 1)
	y = matrix_a[:, 0].reshape(n + 1, n + 1)
	im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
	im.set_rasterized(True)
	return im

def format_ticks(ax, x=True, y=True, rotate=True):
	bp_formatter = EngFormatter('b')
	if y:
		ax.yaxis.set_major_formatter(bp_formatter)
	if x:
		ax.xaxis.set_major_formatter(bp_formatter)
		ax.xaxis.tick_bottom()
	if rotate:
		ax.tick_params(axis='x',rotation=45)

def data_from_clr(clr, region):
	return clr.matrix(balance=True).fetch(region)

def makeregion(chr, start, windows):
	end = start + 90 * windows[0]
	return (chr, start, end)

def plot_contacts_simple(windows, resolution, data, insulation_table, region, norm):
	plt.rcParams['font.size'] = 12

	f, ax = plt.subplots(figsize=(18, 6))
	im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
	ax.set_aspect(0.5)
	ax.set_ylim(0, 10*windows[0])
	format_ticks(ax, rotate=False)
	ax.xaxis.set_visible(False)

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
	plt.colorbar(im, cax=cax)

	insul_region = bioframe.select(insulation_table, region)

	ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)
	ins_ax.set_prop_cycle(plt.cycler("color", plt.cm.plasma(np.linspace(0,1,5))))
	ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
				insul_region['log2_insulation_score_'+str(windows[0])],
				label=f'Window {windows[0]} bp')

	ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);

	format_ticks(ins_ax, y=False, rotate=False)
	ax.set_xlim(region[1], region[2])
	plt.savefig("fig1.png")

	for res in windows[1:]:
		ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), insul_region[f'log2_insulation_score_{res}'], label=f'Window {res} bp')
	ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);
	# f
	plt.savefig("fig2.png")

def plot_with_boundaries(insulation_table, region, data, resolution, norm, windows):
	f, ax = plt.subplots(figsize=(20, 10))
	im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
	ax.set_aspect(0.5)
	ax.set_ylim(0, 10*windows[0])
	format_ticks(ax, rotate=False)
	ax.xaxis.set_visible(False)

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
	plt.colorbar(im, cax=cax)

	insul_region = bioframe.select(insulation_table, region)

	ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)

	ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
				insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

	boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[0]}'])]
	weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']]
	strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']]
	ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
				weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
	ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
				strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

	ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);

	format_ticks(ins_ax, y=False, rotate=False)
	ax.set_xlim(region[1], region[2])
	plt.savefig("chr_with_boundaries.png")

def make_histkwargs():
	histkwargs = dict(
		bins=10**np.linspace(-4,1,200),
		histtype='step',
		lw=2,
	)
	return histkwargs

def boundary_strength_and_plot(insulation_table, windows, histkwargs):
	f, axs = plt.subplots(len(windows),1, sharex=True, figsize=(6,6), constrained_layout=True)
	for i, (w, ax) in enumerate(zip(windows, axs)):
		ax.hist(
			insulation_table[f'boundary_strength_{w}'],
			**histkwargs
		)
		ax.text(0.02, 0.9,
				 f'Window {w//1000}kb',
				 ha='left',
				 va='top',
				 transform=ax.transAxes)

		ax.set(
			xscale='log',
			ylabel='# boundaries'
		)

	axs[-1].set(xlabel='Boundary strength');
	plt.savefig("boundary_strength.png")

def call_boundaries(insulation_table, windows, histkwargs):
	f, axs = plt.subplots(len(windows), 1, sharex=True, figsize=(6,6), constrained_layout=True)
	thresholds_li = {}
	thresholds_otsu = {}
	for i, (w, ax) in enumerate(zip(windows, axs)):
		ax.hist(
			insulation_table[f'boundary_strength_{w}'],
			**histkwargs
		)
		thresholds_li[w] = threshold_li(insulation_table[f'boundary_strength_{w}'].dropna().values)
		thresholds_otsu[w] = threshold_otsu(insulation_table[f'boundary_strength_{w}'].dropna().values)
		n_boundaries_li = (insulation_table[f'boundary_strength_{w}'].dropna()>=thresholds_li[w]).sum()
		n_boundaries_otsu = (insulation_table[f'boundary_strength_{w}'].dropna()>=thresholds_otsu[w]).sum()
		ax.axvline(thresholds_li[w], c='green')
		ax.axvline(thresholds_otsu[w], c='magenta')
		ax.text(0.01, 0.9,
				 f'Window {w//1000}kb',
				 ha='left',
				 va='top',
				 transform=ax.transAxes)
		ax.text(0.01, 0.7,
				f'{n_boundaries_otsu} boundaries (Otsu)',
				c='magenta',
				ha='left',
				va='top',
				transform=ax.transAxes)
		ax.text(0.01, 0.5,
				f'{n_boundaries_li} boundaries (Li)',
				c='green',
				ha='left',
				va='top',
				transform=ax.transAxes)

		ax.set(
			xscale='log',
			ylabel='# boundaries'
		)

	axs[-1].set(xlabel='Boundary strength')
	plt.savefig("boundaries.png")
	return thresholds_li, thresholds_otsu

def get_ctcf_data(data_dir):
	ctcf_fc_file = cooltools.download_data("HFF_CTCF_fc", cache=True, data_dir=data_dir)
	return ctcf_fc_file

def tabulate_boundaries(insulation_table, windows):
	is_boundary = np.any([
			~insulation_table[f'boundary_strength_{w}'].isnull()
			for w in windows],
		axis=0)
	boundaries = insulation_table[is_boundary]
	# boundaries.head()
	return boundaries

def chip_signal(boundaries, data_dir):
	# Calculate the average ChIP singal/input in the 3kb region around the boundary.
	flank = 1000 # Length of flank to one side from the boundary, in basepairs
	ctcf_chip_signal = bbi.stackup(
		data_dir+'/test_CTCF.bigWig',
		boundaries.chrom,
		boundaries.start-flank,
		boundaries.end+flank,
		bins=1).flatten()
	return ctcf_chip_signal

def plot_boundary_strength_vs_ctcf_enrich(boundaries, windows, ctcf_chip_signal, thresholds_li, thresholds_otsu):
	w=windows[0]
	f, ax = plt.subplots()
	ax.loglog(
		boundaries[f'boundary_strength_{w}'],
		ctcf_chip_signal,
		'o',
		markersize=1,
		alpha=0.05
	);
	ax.set(
		xlim=(1e-4,1e1),
		ylim=(3e-2,3e1),
		xlabel='Boundary strength',
		ylabel='CTCF enrichment over input')

	ax.axvline(thresholds_otsu[w], ls='--', color='magenta', label='Otsu threshold')
	ax.axvline(thresholds_li[w], ls='--', color='green', label='Li threshold')
	ax.legend()
	plt.savefig("boundaries_vs_ctcf.png")

def threshold_and_plot_boundaries_with_ctcf(boundaries, ctcf_chip_signal, histkwargs, windows):
	f, ax = plt.subplots()
	ax.set(xscale='log', xlabel='Boundary strength')
	ax.hist(
		boundaries[f'boundary_strength_{windows[0]}'][ctcf_chip_signal>=2],
		label='CTCF Chip/Input ≥ 2.0',
		**histkwargs
	);
	ax.hist(
		boundaries[f'boundary_strength_{windows[0]}'][ctcf_chip_signal<2],
		label='CTCF Chip/Input < 2.0',
		**histkwargs
	);
	ax.hist(
		boundaries[f'boundary_strength_{windows[0]}'],
		label='all boundaries',
		**histkwargs
	);
	ax.legend(loc='upper left')
	plt.savefig("ctcf_boundaries.png")

def pileup_1d_select_boundaries(insulation_table, windows, thresholds_otsu):
	# Select the strict thresholded boundaries for one window size
	top_boundaries = insulation_table[insulation_table[f'boundary_strength_{windows[1]}']>=thresholds_otsu[windows[1]]]
	return top_boundaries

def pileup_1d_create_stackup(data_dir, top_boundaries, resolution):
	# Create of the stackup, the flanks are +- 50 Kb, number of bins is 100 :
	flank = 50000 # Length of flank to one side from the boundary, in basepairs
	nbins = 100   # Number of bins to split the region
	stackup = bbi.stackup(data_dir+'/test_CTCF.bigWig',
						  top_boundaries.chrom,
						  top_boundaries.start+resolution//2-flank,
						  top_boundaries.start+resolution//2+flank,
						  bins=nbins)
	return stackup, flank, nbins

def pileup_1d_plot(stackup, flank, nbins):
	f, ax = plt.subplots(figsize=[7,5])
	ax.plot(np.nanmean(stackup, axis=0) )
	ax.set(xticks=np.arange(0, nbins+1, 10),
		   xticklabels=(np.arange(0, nbins+1, 10)-nbins//2)*flank*2/nbins/1000,
		   xlabel='Distance from boundary, kbp',
		   ylabel='CTCF ChIP-Seq mean fold change over input');
	plt.savefig("pileup.png")

def main():
	print("hi")
	data_dir, resolution, clr, windows, insulation_table, norm = getdata()
	print_first_summary(windows, insulation_table)

	region = makeregion("chr2", 10_500_000, windows)
	data = data_from_clr(clr, region)

	plot_contacts_simple(windows, resolution, data, insulation_table, region, norm)
	plot_with_boundaries(insulation_table, region, data, resolution, norm, windows)

	histkwargs = make_histkwargs()
	boundary_strength_and_plot(insulation_table, windows, histkwargs)
	thresholds_li, thresholds_otsu = call_boundaries(insulation_table, windows, histkwargs)

	ctcf_fc_file = get_ctcf_data(data_dir)

	boundaries = tabulate_boundaries(insulation_table, windows)
	ctcf_chip_signal = chip_signal(boundaries, data_dir)
	plot_boundary_strength_vs_ctcf_enrich(boundaries, windows, ctcf_chip_signal, thresholds_li, thresholds_otsu)
	threshold_and_plot_boundaries_with_ctcf(boundaries, ctcf_chip_signal, histkwargs, windows)

	top_boundaries = pileup_1d_select_boundaries(insulation_table, windows, thresholds_otsu)
	stackup, flank, nbins = pileup_1d_create_stackup(data_dir, top_boundaries, resolution)
	pileup_1d_plot(stackup, flank, nbins)


if __name__ == "__main__":
	main()
