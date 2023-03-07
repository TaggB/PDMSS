#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 20:12:54 2022

Converts .ab1 files to .csv of structure:
    row = base number
    columns = channel1 read,c hannel2 read, channel3 read, channel4 read, called nucleotide, phred score

requirements:
    - pandas. Exists on base anaconda installation


- uses code from http://github.com/bow/abifpy to read in file initially
  as noted in section 2


Structure:
    #section 1: code to convert ab1 to csv & perform downstream processing
    #section2: ab1 reader code (modified by @benjamintagg)
    #section3: @main arg: allows run from command line


@author: benjamintagg
"""
# =============================================================================
# Section 1: author @benjamintagg
# =============================================================================
# config
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import copy
wvdict = {540:'ddGTP'}
wvdict.update({568:'ddATP'})
wvdict.update({595:'ddTTP'})
wvdict.update({615:'ddCTP'})
###

def ab1_to_csv(file):
    # for wavelengths given in wv vars: 
        # a = 568
        # c = 615
        # t= 595
        # g = 540
    handle = Trace(file,trimming=False)
    gt_seq = list(handle.seq) # base called based on phred
    #seq_start = handle.seq[1] # deprecated in below (by request of IC to get complete read)
    #seq_end = handle.seq[2]

    qual_val = list(handle.qual_val) # phred qual value, derived from data below
   
   # using phred values, and raw data below in DATA9:12 inclusive
    # a trace illustrating bases cna be shown by:
        #plt.hist(sequence1,bins=len(sequence1))
        # and then repeating for sequence 2:4

    sequence1 = pd.Series(handle.get_data('DATA9'))#[seq_start:seq_end]) # fluorescences at all time points
    channel1 = wvdict[handle.get_data('DyeW1')]
    sequence2 = pd.Series(handle.get_data('DATA10'))#[seq_start:seq_end])
    channel2 = wvdict[handle.get_data('DyeW2')]
    sequence3 = pd.Series(handle.get_data('DATA11'))#[seq_start:seq_end])
    channel3 = wvdict[handle.get_data('DyeW3')]
    sequence4 = pd.Series(handle.get_data('DATA12'))#[seq_start:seq_end])
    channel4 = wvdict[handle.get_data('DyeW4')]
    
    posn = handle.get_data('PLOC1')#[seq_start:seq_end] # position for base call, trimmed
    
    
    fluors1 = sequence1.iloc[list(posn)] # fluoresences at base call time points
    fluors2 = sequence2.iloc[list(posn)]
    fluors3 = sequence3.iloc[list(posn)]
    fluors4 = sequence4.iloc[list(posn)]
    
    # create arrays for qual score and base calls (which are peak measures), and pad with zeros
    # to be equal in length to data from all time points
    # for base calls, other time  points are padded with "N", but should be distinguished from "N" calls at peaks
    all_seqs = np.chararray(np.size(sequence1))
    all_phred = np.zeros([np.size(sequence1)])
    all_seqs.fill("N")
    for base_call_no, flor_num in enumerate(fluors1.index):
        all_seqs[flor_num] = gt_seq[base_call_no]
        all_phred[flor_num] = qual_val[base_call_no]

    
    # put into pandas DataFrame
    #seqdf = pd.concat((sequence1,sequence2,sequence3,sequence4,gt_seq,qual_val),axis=1)
    # added all time points
    seqdf = pd.concat((sequence1,sequence2,sequence3,sequence4,fluors1,fluors2,fluors3,fluors4,pd.Series(all_phred)),axis=1)
    seqdf['seq'] = all_seqs.decode(encoding="UTF-8")
    
    
    # here: issue is that seq being lost on concat
    # created as series and adds a 'b' - insepct test
    
    
    seqdf.columns = [channel1,channel2,channel3,channel4,channel1+"peak",channel2+"peak",channel3+"peak",channel4+"peak",'phred score','seq']
    savepath = "/".join(file.split('/')[:-1])+"/"+"".join(file.split('/')[-1].split('.')[:-1])+"_raw.csv"
    seqdf.to_csv(savepath)
    return(seqdf,savepath)

def do_beat(seqdf,savepath,window = [20,20]):
    """
    Window: the regions to remove before background base noise calculation. 
        Default is first 100 bases and last 50.
    If trimming was performed with phred, window =[20,20]

    """
    # trim beginning and end
    trimmed_df = seqdf.iloc[window[0]:-window[1]]
    
    # then take values where base != base and calculated median absolute deviation for each base
    # for base in col 0:
    # for base in col0, expected base given by col0_base (and etc for other cols)
    col0_N = trimmed_df[trimmed_df.iloc[:,4]!=trimmed_df.iloc[:,0].name[-3]]
    col1_N = trimmed_df[trimmed_df.iloc[:,4]!=trimmed_df.iloc[:,1].name[-3]]
    col2_N = trimmed_df[trimmed_df.iloc[:,4]!=trimmed_df.iloc[:,2].name[-3]]
    col3_N = trimmed_df[trimmed_df.iloc[:,4]!=trimmed_df.iloc[:,3].name[-3]]

    # first, remove outliers:
    # Get MAD: get absolute (median - xi), and take median = MAD
    # where base call != base in column x
    col0_N_MAD = ((col0_N.iloc[:,0] - col0_N.iloc[:,0].median(axis=0)).abs()).median()
    col1_N_MAD = ((col1_N.iloc[:,1] - col1_N.iloc[:,1].median(axis=0)).abs()).median()
    col2_N_MAD = ((col2_N.iloc[:,2] - col2_N.iloc[:,2].median(axis=0)).abs()).median()
    col3_N_MAD = ((col3_N.iloc[:,3] - col3_N.iloc[:,3].median(axis=0)).abs()).median()
    # when lots of zeros, median = 0
    
    # use MAD to remove outliers (standard practice is to take values within 2.5*MAD +- median)
    if col0_N_MAD > 0:
        if col0_N.iloc[:,0].median(axis=0) >0:
            col0 = trimmed_df.iloc[:,0][(trimmed_df.iloc[:,0]<=((trimmed_df.iloc[:,0].median(axis=0))+(2.5*col0_N_MAD))) & (trimmed_df.iloc[:,0]>=((trimmed_df.iloc[:,0].median(axis=0))-(2.5*col0_N_MAD)))]
        else:
            col0 = trimmed_df.iloc[:,0]
    else:
        col0 = trimmed_df.iloc[:,0]
    if col1_N_MAD !=0:
        if col1_N.iloc[:,0].median(axis=0) !=0:
            col1 = trimmed_df.iloc[:,1][(trimmed_df.iloc[:,1]<=((trimmed_df.iloc[:,1].median(axis=0))+(2.5*col1_N_MAD))) & (trimmed_df.iloc[:,1]>=((trimmed_df.iloc[:,1].median(axis=0))-(2.5*col1_N_MAD)))]
        else:
            col1 = trimmed_df.iloc[:,1]
    else:
        col1 = trimmed_df.iloc[:,1]
    if col2_N_MAD !=0:
        if col2_N.iloc[:,0].median(axis=0) !=0:
            col2 = trimmed_df.iloc[:,2][(trimmed_df.iloc[:,2]<=((trimmed_df.iloc[:,2].median(axis=0))+(2.5*col2_N_MAD))) & (trimmed_df.iloc[:,2]>=((trimmed_df.iloc[:,2].median(axis=0))-(2.5*col2_N_MAD)))]
        else:
            col2 = trimmed_df.iloc[:,2]
    else:
        col2 = trimmed_df.iloc[:,2]
    if col3_N_MAD !=0:
        if col3_N.iloc[:,0].median(axis=0) !=0:
            col3 = trimmed_df.iloc[:,3][(trimmed_df.iloc[:,3]<=((trimmed_df.iloc[:,3].median(axis=0))+(2.5*col3_N_MAD))) & (trimmed_df.iloc[:,3]>=((trimmed_df.iloc[:,3].median(axis=0))-(2.5*col3_N_MAD)))]
        else:
            col3 = trimmed_df.iloc[:,3]
    else:
        col3 = trimmed_df.iloc[:,3]

    # these should share an index, so cna concatenate
    trimmed_df = pd.concat((col0,col1,col2,col3),axis=1)
    # issue above is when MAD = 0, col = 0 (hence added catch for =0)
    # and whne median = 0, col = 0 (added catch for that)
    
    # then perform backgorudn subtraction
        # in beat, they take average as background nosie for each base
        # and subtract it.
        # since read can be 0, this is what offsets peaks, e.g.:
    trimmed_df = trimmed_df -(trimmed_df.mean(axis=0))
        # can be plotted after import matplotlib.pyplot as plt.
    #plt.plot(background_subtracted.index,test.iloc[:,0])
    #plt.plot(background_subtracted.index,test.iloc[:,1])
    #plt.plot(background_subtracted.index,test.iloc[:,2])
    #plt.plot(background_subtracted.index,test.iloc[:,3])
    
    # this is a shit approach - better ways to remove outliers.
    # make <0 (i.e. backgorund = 0)
    trimmed_df[trimmed_df<=0]=0
    # then calculate percentage of each base at each position
    df_pct = pd.concat((trimmed_df.iloc[:,0]/trimmed_df.sum(axis=1),trimmed_df.iloc[:,1]/trimmed_df.sum(axis=1),trimmed_df.iloc[:,2]/trimmed_df.sum(axis=1),trimmed_df.iloc[:,3]/trimmed_df.sum(axis=1)),axis=1)
    df_pct.columns = seqdf.columns[:4]
    # and then can concat base name and phred score
    BEAT_return = pd.concat((df_pct,seqdf.iloc[:,-2:]),axis=1)
    BEAT_return.to_csv(savepath+"BEAT.csv")
    return()


def peak_solver_sorter(seqdf,savepath,threshold = 200,timewindow = 3, trimwindow = [20,20], background = 10):
    """
    Identifies any peaks not used in base call
        Secondary peaks are identified as having fluorescence >= 40 (based on distribution of fluorescence values)
        But independently on each channel
        If they are within +- timewindow of an existing peak, then a max operation is applied over the window
        If they are not, a new peak is added
    Performs local background subtraction
    And quantifies the ratio of bases at each peak.
        N.b. ratios cna be 
        
    A saved PISS file contains backgorund subtracted fluorescence intensities for all peaks,
    the fraction of DNA for each base at that point
    As well as the base call based on threshold.
    
    Runtime is proportional to size of timewindow and background
        
    
    threhsold: the intensity threshold. If two bases have intesnity > threshold, the base at that call will be "N"
    timewindow: the window over which closest peaks are combined (+timewindow,-timewindow) by taking their maximum
    trimwindow: the first x and alst y of [x,y] amino acids that are removed before beat analysis is performed
    background: the length of residues over which an average is performed to determine background for each channel across those residues

    """
    # identify all peaks
    ch1peaks = find_peaks(seq_df.iloc[:,0],height=40)[0]
    ch2peaks = find_peaks(seq_df.iloc[:,1],height=40)[0]
    ch3peaks = find_peaks(seq_df.iloc[:,2],height=40)[0]
    ch4peaks = find_peaks(seq_df.iloc[:,3],height=40)[0]
    all_peaks = np.concatenate((ch1peaks,ch2peaks,ch3peaks,ch4peaks))
    # get locations already identified as peaks:
    existing_calls = seq_df.iloc[:,5].index[~pd.isna(seq_df.iloc[:,5])]
    
    
    ##### had just changed the above line and now won't run
    
    # if not already identified as a peak,  add to an variable  new peaks
    new_peaks = pd.DataFrame()
    for item in all_peaks:
        # if item in existing_calls is fine
        if item not in existing_calls:
            new_peaks = pd.concat((new_peaks,seq_df.iloc[item,:4]),axis=1)
    new_peaks = new_peaks.transpose()
    
    # go through new peaks, find closest existing peaks
    # if there is an existing peak within (timewindow) timepoints of then newpeak (set by timewindow)
        # take a maximum over each channel for newpeak and existing peak
    # otherwise:
        # add a new peak
    combined_peaks = pd.DataFrame()
    original_peaks = []
    for item in new_peaks.index:
        loc_diff = existing_calls - item
        loc_diff.index = np.arange(np.size(loc_diff)) # this is sloppy and inefficient
        # find if any close peaks
        if np.any((loc_diff<=timewindow) & (loc_diff>=-timewindow)):
            # if there are, identify them (can be multiple, so take closest)
           # all_local_peaks = loc_diff.index[(loc_diff<=timewindow) & (loc_diff>=-timewindow)]
            closest_peak = existing_calls[loc_diff.index[loc_diff.index[(loc_diff<=timewindow) & (loc_diff>=-timewindow)]][loc_diff[loc_diff.index[(loc_diff<=timewindow) & (loc_diff>=-timewindow)]]==loc_diff[loc_diff.index[(loc_diff<=timewindow) & (loc_diff>=-timewindow)]].min()]]
            # append this to the list original peaks to retain time indexing (added back later)
            original_peaks.append(closest_peak[0])
            # and get max vlaue across both old and new peak
            try: # lazy fix for when both conditions met
                combined_peak = (pd.concat((seq_df.loc[closest_peak].iloc[:,:4],pd.DataFrame(new_peaks.loc[item]).transpose()),axis=0)).max(axis=0)
            except pd.errors.InvalidIndexError:
                combined_peak = (pd.concat((seq_df.loc[closest_peak].iloc[:,:4],pd.DataFrame(new_peaks.loc[item].iloc[0,:]).transpose()),axis=0)).max(axis=0)
            # and set position = location original peak (done by index later)
        else:
            # if not in existing calls, add new peak
            combined_peak = new_peaks.loc[item]
                # and append its location
            original_peaks.append(item)
        # append to the dataframe containing all new peaks and combined local peaks
        combined_peaks = pd.concat((combined_peaks,combined_peak),axis=1)
        

    combined_peaks = combined_peaks.transpose()
    combined_peaks.index = original_peaks
        
    # then need to get any exsiting peaks not in combined peaks
    for item in existing_calls:
        if item not in combined_peaks.index:
            combined_peaks = pd.concat((combined_peaks,seqdf.loc[existing_calls].iloc[:,:4]),axis=0)
            
    # reorder so that position is incremental

    combined_peaks.sort_index(inplace=True)
    # de-duplicate
    combined_peaks = combined_peaks[~combined_peaks.index.duplicated(keep='first')]
    # rename cols
    combined_peaks.columns = seq_df.iloc[:,:4].columns
    
    # call a base for each position
    # using max
    # if more than two bases > 0, call N
    base_calls = []
    for item in combined_peaks.index:
        if ((combined_peaks.loc[item]>threshold).sum()) >1:
            base_calls.append("N")
        else:
            base_calls.append(combined_peaks.columns[combined_peaks.loc[item]==combined_peaks.loc[item].max()][0][-3])
    #add base calls
    combined_peaks['seq'] = base_calls # combined peaks contains only peak info,
    # because of the max over timewindow
    
    # then perform beat, or don't because it's shit
    # beat_off(seqdf, savepath) # this won't currently run since seqdf is not as expected

    
    # but basically, for backgorund subtraction, couple of thoughts:
        # 1. really important to get samples right - probably want a few different primers
            # and independent sequencing reactions of the same sample, at least
        # 2. have rolling average over bases for background
        # 2. (alternative): use avg of where a base is not the base call as background
        # 2 (alternative): for flip/flop, geenerate background using bases they have in common
            # and apply to the bases where they are different.
            
    # first, trim bases
    combined_peaks = combined_peaks.loc[combined_peaks.index[(combined_peaks.index>=trimwindow[0]) & (combined_peaks.index<=combined_peaks.index.max()-trimwindow[1])]]
    
    #trying strategy of using where base x is not the called base to determine background
    # this seems to work quite nicely when background variable = 10
    
    # laziest way to do this is to go bidirectionally, and then average across two vectors for background
    # even lazier way is to assume that last few bases (i.e. size(seq) - (backrepeats * bacgkround) are not of interest
    
    # for each channel, find instances where base x is not the called base
    # excluding instances where N is called (i.e. where could be multiple peaks)
    bckch0 = combined_peaks.iloc[:,0][(combined_peaks.iloc[:,4]!=combined_peaks.columns[0][2]) & (combined_peaks.iloc[:,4]!='N')]
    bckch1 = combined_peaks.iloc[:,1][(combined_peaks.iloc[:,4]!=combined_peaks.columns[1][2]) & (combined_peaks.iloc[:,4]!='N')]
    bckch2 = combined_peaks.iloc[:,2][(combined_peaks.iloc[:,4]!=combined_peaks.columns[2][2]) & (combined_peaks.iloc[:,4]!='N')]
    bckch3 = combined_peaks.iloc[:,3][(combined_peaks.iloc[:,4]!=combined_peaks.columns[3][2]) & (combined_peaks.iloc[:,4]!='N')]
    allbck = pd.concat((bckch0,bckch1,bckch2,bckch3),axis=1)

    # basecalls are retained as were
    background_subtracted = pd.DataFrame()
    back_repeats = np.floor(combined_peaks.index.max()/background)
    for repeat in np.arange(0,back_repeats-1):
        if repeat == 0:
            peakrange = np.arange(background).astype(int)
        else:
            peakrange = np.arange((repeat)*background,(repeat+1)*background).astype(int)
        #if allbck.index.intersection(peakrange).size !=0: # lazy catch. repeated below
        # if there is backgorund in this range
        # if there are peaks in this range
        if allbck.loc[allbck.index.intersection(peakrange)].size !=0:
            # if there is background, subtract it
            window_avg_back = (allbck.loc[allbck.index.intersection(peakrange)].mean(axis=0)).transpose()
            window_avg_back[pd.isna(window_avg_back)] = 0
            subtr = pd.concat((combined_peaks.loc[combined_peaks.index.intersection(peakrange)].iloc[:,:4] - window_avg_back.transpose(),combined_peaks.loc[combined_peaks.index.intersection(peakrange)].iloc[:,4]),axis=1)
            background_subtracted = pd.concat((background_subtracted,subtr),axis=0)
        else: # no background
            subtr = pd.concat((combined_peaks.loc[combined_peaks.index.intersection(peakrange)].iloc[:,:4],combined_peaks.loc[combined_peaks.index.intersection(peakrange)].iloc[:,4]),axis=1)
            background_subtracted = pd.concat((background_subtracted,subtr),axis=0)
        
    # probably worth noting when background is subspiciously high - shouldn't ever be the case though.
    
    # for final steps of beat equivalent (to get ratios),
    # supposed to remove outliers. Not currently done
    # then get ratios
    forbeat = copy.copy(background_subtracted.iloc[:,:4])
    # zero the background 
    forbeat[forbeat<0] = 0
    forbeat[forbeat==0] = np.nan # for ease
    axisdiv = lambda x: x/x.sum()
    forbeat = forbeat.apply(axisdiv,axis=1)
    forbeat[pd.isna(forbeat)] = 0

    returns = pd.concat((background_subtracted.iloc[:,:4],forbeat,background_subtracted.iloc[:,4]),axis=1)
    returns.to_csv(savepath+"PISS.csv")

    # here: nb that bits after sequence cna have negativ ebeat
    
    return()
    


# =============================================================================
#  section 2: one line was modified to make data comparable to seq (annotate)
    # otherwise, code was taken from


#Abfipy reader file 
#: from http://github.com/bow/abifpy
# circumvents biopython requirement
# =============================================================================
import datetime
import struct
from os.path import splitext, basename

from sys import version_info

RELEASE = False
__version_info__ = ('1', '0', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''


__all__ = ['Trace']

# dictionary for deciding which values to extract and contain in self.data
EXTRACT = {
            'TUBE1': 'well',
            'DySN1': 'dye',
            'GTyp1': 'polymer',
            'MODL1': 'model', 
            'RUND1': 'run start date',
            'RUND2': 'run finish date',
            'RUND3': 'data collection start date',
            'RUND4': 'data collection finish date',
            'RUNT1': 'run start time',
            'RUNT2': 'run finish time',
            'RUNT3': 'data collection start time',
            'RUNT4': 'data collection finish time',
            'DATA1': 'raw1',
            'DATA2': 'raw2',
            'DATA3': 'raw3',
            'DATA4': 'raw4',
            'PLOC2': 'tracepeaks',
            'FWO_1': 'baseorder',
          }     

# dictionary for unpacking tag values
_BYTEFMT = {
            1: 'b',     # byte
            2: 's',     # char
            3: 'H',     # word
            4: 'h',     # short
            5: 'i',     # long
            6: '2i',    # rational, legacy unsupported
            7: 'f',     # float
            8: 'd',     # double
            10: 'h2B',  # date
            11: '4B',   # time
            12: '2i2b', # thumb
            13: 'B',    # bool
            14: '2h',   # point, legacy unsupported
            15: '4h',   # rect, legacy unsupported
            16: '2i',   # vPoint, legacy unsupported
            17: '4i',   # vRect, legacy unsupported
            18: 's',    # pString
            19: 's',    # cString
            20: '2i',   # Tag, legacy unsupported
           }

# header structure
_HEADFMT = '>4sH4sI2H3I'

# directory data structure
_DIRFMT = '>4sI2H4I'

# to handle py3 IO
def py3_get_string(byte):
    if version_info[0] < 3:
        return byte
    else:
        return byte.decode()

def py3_get_byte(string):
    if version_info[0] < 3:
        return string
    else:
        return string.encode()

class Trace(object):
    """Class representing trace file."""
    def __init__(self, in_file, trimming=False):        
        self._handle = open(in_file, 'rb')
        try:
            self._handle.seek(0)
            if not self._handle.read(4) == py3_get_byte('ABIF'):
                raise IOError('Input is not a valid trace file')
        except IOError:
            self._handle = None
            raise
        else:
            # header data structure:
            # file type, file, version, tag name, tag number, element type code,
            # element size, number of elements, data size, data offset, handle,
            # file type, file version
            # dictionary for containing file metadata
            self.data = {}
            # dictionary for containing extracted directory data
            self.tags = {}
            self.trimming = trimming
            # values contained in file header
            self._handle.seek(0)
            header = struct.unpack(_HEADFMT, 
                     self._handle.read(struct.calcsize(_HEADFMT)))
            # file format version
            self.version = header[1]

            # build dictionary of data tags and metadata
            for entry in self._parse_header(header):
                key = entry.tag_name + str(entry.tag_num)
                self.tags[key] = entry
                # only extract data from tags we care about
                if key in EXTRACT:
                    # e.g. self.data['well'] = 'B6'
                    self.data[EXTRACT[key]] = self.get_data(key)

            self.id = self._get_file_id(in_file)
            self.name = self.get_data('SMPL1')
            self.seq = self.get_data('PBAS2')
            self.qual = ''.join([chr(ord(value) + 33) for value in self.get_data('PCON2')])
            self.qual_val = [ord(value) for value in self.get_data('PCON2')]

            if trimming:
                self.seq, self.qual, self.qual_val = map(self.trim, 
                                                        [self.seq, self.qual,
                                                        self.qual_val])

    def __repr__(self):
        """Represents data associated with the file."""
        if len(self.seq) > 10:
            seq = "{0}...{1}".format(self.seq[:5], self.seq[-5:])
            qual_val = "[{0}, ..., {1}]".format(
                      repr(self.qual_val[:5])[1:-1], 
                      repr(self.qual_val[-5:])[1:-1])
        else:
            seq = self.seq
            qual_val = self.qual_val

        return "{0}({1}, qual_val:{2}, id:{3}, name:{4})".format(
                self.__class__.__name__, repr(seq), qual_val,
                repr(self.id), repr(self.name))
    
    def _parse_header(self, header):
        """Generator for directory contents."""
        # header structure:
        # file signature, file version, tag name, tag number, 
        # element type code, element size, number of elements
        # data size, data offset, handle
        head_elem_size = header[5]
        head_elem_num = header[6]
        head_offset = header[8]
        index = 0
        
        while index < head_elem_num:
            start = head_offset + index * head_elem_size
            # added directory offset to tuple
            # to handle directories with data size <= 4 bytes
            self._handle.seek(start)
            dir_entry =  struct.unpack(_DIRFMT, 
                        self._handle.read(struct.calcsize(_DIRFMT))) + (start,)
            index += 1
            yield _TraceDir(dir_entry, self._handle)

    def _get_file_id(self, in_file):
        """Returns filename without extension."""
        return splitext(basename(in_file))[0]

    def close(self):
        """Closes the Trace file object."""
        self._handle.close()
    

    def get_data(self, key):
        """Returns data stored in a tag."""
        return self.tags[key].tag_data

    def seq_remove_ambig(self, seq):
        """Replaces extra ambiguous bases with 'N'."""
        import re
        seq = self.seq
        return re.sub("K|Y|W|M|R|S", 'N', seq)

    def export(self, out_file="", fmt='fasta'):
        """Writes the trace file sequence to a fasta file.
        
        Keyword argument:
        out_file -- output file name (detault 'tracefile'.fa)
        fmt -- 'fasta': write fasta file, 'qual': write qual file, 'fastq': write fastq file
        """
        if out_file == "":
            file_name = self.id
            if fmt == 'fasta':
                file_name += '.fa'
            elif fmt == 'qual':
                file_name += '.qual'
            elif fmt == 'fastq':
                file_name += '.fq'
            else:
                raise ValueError('Invalid file format: {0}.'.format(fmt))
        else:
            file_name = out_file
        
        if fmt == 'fasta':
            contents = '>{0} {1}\n{2}\n'.format(
                        self.id, 
                        self.name, 
                        self.seq)
        elif fmt == 'qual':
            contents = '>{0} {1}\n{2}\n'.format(
                        self.id, 
                        self.name, 
                        ' '.join(map(str, self.qual_val)))
        elif fmt == 'fastq':
            contents = '@{0} {1}\n{2}\n+{0} {1}\n{3}\n'.format(
                        self.id, 
                        self.name, 
                        self.seq, ''.join(self.qual))

        with open(file_name, 'w') as out_file:
            out_file.writelines(contents)

    def trim(self, seq, cutoff=0.05):
        """Trims the sequence using Richard Mott's modified trimming algorithm.
        
        Keyword argument:
        seq -- sequence to be trimmed
        cutoff -- probability cutoff value
        Trimmed bases are determined from their segment score, ultimately
        determined from each base's quality values. 
        
        More on:
        http://www.phrap.org/phredphrap/phred.html
        http://www.clcbio.com/manual/genomics/Quality_trimming.html
        """
        # set flag for trimming
        start = False
        # set minimum segment size
        segment = 20
        trim_start = 0
        
        if len(seq) <= segment:
            raise ValueError('Sequence can not be trimmed because \
                             it is shorter than the trim segment size')
        else:
            # calculate probability back from formula used
            # to calculate phred qual values
            score_list = [cutoff - (10 ** (qual/-10.0)) for 
                         qual in self.qual_val]

            # calculate cummulative score_list
            # if cummulative value < 0, set to 0
            # first value is set to 0 (assumption: trim_start is always > 0)
            running_sum = [0]
            for i in range(1, len(score_list)):
                num = running_sum[-1] + score_list[i]
                if num < 0:
                    running_sum.append(0)
                else:
                    running_sum.append(num)
                    if not start:
                        # trim_start = value when cummulative starts to be > 0
                        trim_start = i
                        start = True

            # trim_finish = index of the highest cummulative value,
            # marking the segment with the highest cummulative score 
            trim_finish = running_sum.index(max(running_sum)) 

            return (seq[trim_start:trim_finish],trim_start,trim_finish) # modified to give trim_start,finish

class _TraceDir(object):
    """Class representing directory content."""
    def __init__(self, tag_entry, handle):
        self.tag_name = py3_get_string(tag_entry[0])
        self.tag_num = tag_entry[1]
        self.elem_code = tag_entry[2]
        self.elem_size = tag_entry[3]
        self.elem_num = tag_entry[4]
        self.data_size = tag_entry[5]
        self.data_offset = tag_entry[6]
        self.data_handle = tag_entry[7]
        self.tag_offset = tag_entry[8]

        # if data size is <= 4 bytes, data is stored inside the directory
        # so offset needs to be changed
        if self.data_size <= 4:
            self.data_offset = self.tag_offset + 20

        self.tag_data = self._unpack(handle)

    def __repr__(self):
        """Represents data associated with a tag."""
        summary = ['tag_name: {0}'.format(repr(self.tag_name))]
        summary.append('tag_number: {0}'.format(repr(self.tag_num)))
        summary.append('elem_code: {0}'.format(repr(self.elem_code)))
        summary.append('elem_size: {0}'.format(repr(self.elem_size)))
        summary.append('elem_num: {0}'.format(repr(self.elem_num)))
        summary.append('data_size: {0}'.format(repr(self.data_size)))
        summary.append('data_offset: {0}'.format(repr(self.data_offset)))
        summary.append('data_handle: {0}'.format(repr(self.data_handle)))
        summary.append('tag_offset: {0}'.format(repr(self.tag_offset)))
        summary.append('tag_data: {0}'.format(repr(self.tag_data)))
       
        return '\n'.join(summary)

    def _unpack(self, handle):
        """Returns tag data"""
        if self.elem_code in _BYTEFMT:
            
            # because ">1s" unpacks differently from ">s"
            num = '' if self.elem_num == 1 else str(self.elem_num)
            fmt = "{0}{1}{2}".format('>', num, _BYTEFMT[self.elem_code])
            start = self.data_offset
    
            handle.seek(start)
            data = struct.unpack(fmt, handle.read(struct.calcsize(fmt)))
            
            # no need to use tuple if len(data) == 1
            if self.elem_code not in [10, 11] and len(data) == 1:
                data = data[0]

            # account for different data types
            if self.elem_code == 2:
                return py3_get_string(data)
            elif self.elem_code == 10:
                return datetime.date(*data)
            elif self.elem_code == 11:
                return datetime.time(*data)
            elif self.elem_code == 13:
                return bool(data)
            elif self.elem_code == 18:
                return py3_get_string(data[1:])
            elif self.elem_code == 19:
                return py3_get_string(data[:-1])
            else:
                return data
        else:
            return None
if __name__ == "__main__":
    file = input("Drag in an ab1 file, or enter its filepath. Then press enter >")
    seq_df,savepath = ab1_to_csv(file.strip("'"))
    #do_beat(seq_df,savepath)
    peak_solver_sorter(seq_df,savepath)
    print("\n \n Raw file has been saved to {} as {}".format("/".join(file.split('/')[:-1]),"".join(file.split('/')[-1].split('.')[:-1])+".csv"))
   # print("\n PISS File has been saved to {} as {}".format("/".join(file.split('/')[:-1]),"".join(file.split('/')[-1].split('.')[:-1])+"_PISS.csv"))

    # savepath
    
