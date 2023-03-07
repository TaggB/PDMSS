# PDMSS
## Peak discrimination for mixed sanger sequences
</br>
</br>
The aim of this software is to allow calculation of the relative abundance of two or more DNA sequences from singular Sanger sequencing reads over 
</br>
</br>
Following expression of multiple splice variants of a gene, throught RT-PCR, we might obtain snager sequences that are aligned at some point, but differ in the identity of bases at a given position. We want to quantify the expression of these different genes. This is not trivial due to a few factors:  
</br>
</br>
1. The resolution of peaks changes along the sequencing read. Initially, it is very poor, since sequencing primers are often not fully bound. The same occurs at the end of the read. We can trim these regions.  
2. Each of four florescence channels (ddATP, ddTTP, ddCTP, ddGTP) has its own 'background' fluorescence level, on top of which, the intensity of the signal varies.  
3. C and G regions can cause emergence of 'lag' in the Sanger read, such that peaks are offset (in time) slightly following C/G-rich regions.  
4. Where the bases are different between two channels, we need to identify (amongst all of the above factors), the real difference in their intensity that is attributable to differences in expression of two alternatively spliced variants.  
</br>
</br>
One approach has been proposed for this during short (20 bases) sequence stretches (https://hanlab.cc/beat/), but sequencing reads can be ~1000 base pairs long, for which this procedure does not perform stably. A modified version of that method for long DNA bases is included in the file under the function "do_beat()", but is not called by the __main__ argument.
</br>
</br>
PDMSS offers an alternative. It works in the following way:  
1. All peaks are identified that have amplitude greater than a given threshold (default 40 based on distributions during validation).  
2. These peaks, which may represent bona fide (additional peaks) are merged with those called during Sanger sequencing.  
3. Peaks are combined if there is another peak within a given time window. Take the maximum for that channel during that window to create a new singular peak.  
4. ***Base calls are amde at this point (presence of multiple channels = N)***. Because of this, base calls might not match the expectation.
5. Trim the sequence to remove low quality read sequences.  
6. Perform background subtraction for each channel separately, using a rolling window to account for changing quality of read.  
7. Normalise each channel to max amplitude.  
8. Calculate the fraction of signal at that position attributable to a given base. For a perfect call, a channel would have a score of 1.0, and all other channels 0.0. A mixture of DNA will have values between this.
</br>
</br>
## Usage
</br>
This software uses some code from http://github.com/bow/abifpy to load .abi files, in accordance with fair use.
</br>
Otherwise, this can be run from terminal using a stable version of anaconda, and can be called from terminal by running the file and dragging in a .abi file when prompted. Two files are returned: one containing the ***full** raw sanger sequence as .csv, and another containing the PDMSS sequence as a .csv.  
</br>
There are a few parameters to play with that should be validated (and can be modified in the kwargs).
These are idenified in the function description. To read, call help(peak_solver_sorter).
</br>
Of course, plotting and further downstream analysis cna be performed using the .csv file in Python using standard libraries, such as matplotlib etc. Loading cna be performed using pandas. IPython is recommended for that.
