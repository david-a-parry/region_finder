# region_finder

A python library for parsing and searching genomic intervals in memory. Also provides a means for randomly sampling intervals.

## Installation

Installation using pip:

    python3 -m pip install git+git://github.com/david-a-parry/region_finder.git --user

Or you can clone this repository and install:

    git clone https://github.com/david-a-parry/region_finder.git
    cd region_finder
    python3 setup.py install --user

## Synopsis

### Searching Intervals

Parse a BED file and read into memory. Then create a region finder object to search intervals.

```
>>> from region_finder.bed_parser import BedParser
>>> from region_finder.region_finder import RegionFinder

>>> bed_intervals = BedParser("test/test_data/test_bed.gz")
>>> bed_searcher = RegionFinder(bed_intervals)
```

You can search either using a chromosome, start and end arguments or using an interval string. The two examples below are equivalant (coordinates for searching are 1-based):

```
>>> results = bed_searcher.fetch("20", 674880, 674916)
>>> results = bed_searcher.fetch_by_interval("20:674880-674916")
```

Results are lists of GenomicInterval objects:

```
>>> results
[<region_finder.genomic_interval.GenomicInterval object at 0x7fd0d8f2cfc0>, <region_finder.genomic_interval.GenomicInterval object at 0x7fd0d8f2d080>]
>>> results[0].contig
'20'
>>> results[0].start  # start coordinates are 0-based like BED regions
674693
>>> results[0].end  # end coordinates are 1-based like BED regions
674883
>>> str(results[0])
'20:674694-674883'
>>> str(results[1])
'20:674884-675056'
```

GenomicInterval objects may contain multiple overlapping intervals from original file. Access the original intervals using the 'regions' attribute:

```
>>> results[1].regions
[['20', 674883, 674916, 'MLT1E3', '759', '-'], ['20', 674915, 675056, 'MLT1E3', '451', '-']]
```

For more information about Genomic Interval objects:

```
>>> from region_finder.genomic_interval import GenomicInterval
>>> help(GenomicInterval)
```

You can pass a list of intervals in standard "\<chr\>:\<start\>-\<end\>" format instead of a BED file using the RegionIter class:

```
>>> from region_finder.region_iter import RegionIter
>>> input_regions = ['20:674383-674693', '20:674694-674883', '20:674884-675056']
>>> reg_iter = RegionIter(input_regions)
>>> reg_searcher = RegionFinder(reg_iter)
>>> results = reg_searcher.fetch('20', 674700, 674800)
>>> str(results[0])
'20:674694-674883'
>>> results[0].regions
[['20', 674693, 674883]]
```

Or you can use the IntervalIter class to process your intervals - an iterable of lists should be provided where each list simply requires that the first three columns correspond to the contig, start (0-based) and end (1-based) of your regions:

```
>>> from region_finder.region_finder import RegionFinder
>>> from region_finder.interval_iter import IntervalIter
>>> intervals = [['20', 674883, 674916, 'MLT1E3', '759', '-'], 
                 ['20', 674915, 675056, 'MLT1E3', '451', '-'],
                 ['20', 674915, 675056, 'Foo', '999', '+']]
>>> intvl_iter = IntervalIter(intervals)
>>> reg_finder = RegionFinder(intvl_iter)
```

### Randomly Sampling Intervals

This module also provides a means for randomly sampling from a set of intervals. Regions are merged and a linear index is created in memory so that any given position is equally likely to be sampled irrespective of whether positions lie within long or short regions or whether positions occur multiple times in overlapping intervals.

As for region searching, the first step is to read your intervals into a BedParser/RegionIter/IntervalIter object where overlapping regions will be merged for sampling:

```
>>> from region_finder.bed_parser import BedParser
>>> from region_finder.interval_sampler import IntervalSampler
>>> bed_intervals = BedParser("test/test_data/test_bed.gz")
>>> sampler = IntervalSampler(bed_intervals)
```

Random (1-based) positions can be sampled:

```
>>> sampler.random_position()
('22', 34976445)
```

Random intervals of a given length can also be sampled:

```
>>> print(sampler.random_interval(10))
20:10286526-10286535
```

You can also ask for a number of randomly sampled intervals using the `random_sample` method. You must provide the mean length and standard deviation of length of regions to sample. The mean and standard deviation will be used to generate a random set of lengths drawn from a normal distribution to return:

```
>>> rnd_regs = sampler.random_sample(n=1000, mean_length=100, sd=50, allow_overlaps=False)
```

The regions returned will be sorted in coordinate order (this facilitates ensuring no overlaps when overlaps are not desired).

Note that in some instances the lengths may deviate slightly from the expected distribution as intervals are truncated or shifted to ensure they lie fully within the original set of intervals. In the above instance where many of the original regions are shorter than 100 bp our mean is lower than expected:

```
>>> import numpy as np
>>> reg_lengths = np.array([len(x) for x in rnd_regs])
>>> reg_lengths.mean()
95.714
>>> reg_lengths.std()
49.200937023597426
```

This is unlikely to be an issue if using a list of much longer regions, but something to be aware of when your requested lengths may exceed the lengths of your original regions.

To retrieve a set of intervals of uniform length, simply set the standard deviation to 0 (caveat from above applies):

```
>>> rnd_regs = sampler.random_sample(n=1000, mean_length=10, sd=0, allow_overlaps=False)
>>> reg_lengths = np.array([len(x) for x in rnd_regs])
>>> reg_lengths.mean()
10.0
>>> reg_lengths.std()
0.0
```

If you know the lengths of regions you want to retrieve, use the random_sample_given_lengths function:

```
>>> my_lengths = [5, 100, 200, 20]
>>> regs = sampler.random_sample_given_lengths(my_lengths)
>>> [len(x) for x in regs]
[5, 200, 20, 100]
>>> # order differs due to coordinate sorting of regions, but each requested length was retrieved
>>> np.sort([len(x) for x in regs]) == np.sort(my_lengths)
array([ True,  True,  True,  True])
```
