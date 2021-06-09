# region_finder

A python library for parsing and searching genomic intervals in memory.

## Synopsis

Parse a BED file and read into memory then create a region finder object to search intervals. 
~~~
>>> from region_finder.bed_parser import BedParser
>>> from region_finder.region_finder import RegionFinder

>>> bed_intervals = BedParser("test/test_data/test_bed.gz")
>>> bed_searcher = RegionFinder(bed_intervals)

~~~

You can search either using a chromosome, start and end arguments or using an interval string. The two examples below are equivalant.
~~~
>>> results = bed_searcher.fetch("20", 674880, 674916)
>>> results = bed_searcher.fetch_by_interval("20:674880-674916")
~~~

Results are lists of GenomicInterval objects

~~~
>>> results
[<region_finder.genomic_interval.GenomicInterval object at 0x7fd0d8f2cfc0>, <region_finder.genomic_interval.GenomicInterval object at 0x7fd0d8f2d080>]
>>> str(results[0])
'20:674694-674883'
>>> str(results[1])
'20:674884-675056'
~~~

GenomicInterval objects may contain multiple overlapping intervals from original file. Access the original intervals using the 'regions' attribute.

~~~
>>> results[1].regions
[['20', 674883, 674916, 'MLT1E3', '759', '-'], ['20', 674915, 675056, 'MLT1E3', '451', '-']]

~~~

For more information about Genomic Interval objects:
~~~
>>> from region_finder.genomic_interval import GenomicInterval
>>> help(GenomicInterval)

~~~

You can pass a list of intervals in standard "\<chr\>:\<start\>-\<end\>" format instead of a BED file.
~~~
>>> input_regions = ['20:674383-674693', '20:674694-674883', '20:674884-675056']
>>> reg_iter = RegionIter(input_regions)
>>> reg_searcher = RegionFinder(reg_iter)
>>> results = reg_searcher.fetch('20', 674700, 674800)
>>> str(results[0])
'20:674694-674883'
>>> results[0].regions
[['20', 674693, 674883]]

~~~
