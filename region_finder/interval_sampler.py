import warnings
import bisect
import numpy as np
from random import randint
from .genomic_interval import GenomicInterval


class IntervalSampler(object):
    '''
    From an IntervalIter object create a linear index of regions and
    facilitate random sampling corrected for region lengths.
    '''
    def __init__(self, interval_iter):
        self.idx = np.cumsum([len(x) for x in interval_iter.intervals],
                             dtype=np.int64)
        self.length = self.idx[-1]
        self.intervals = interval_iter.intervals

    def __len__(self):
        return self.length

    def interval_by_index(self, i):
        if i >= self.length:
            return None
        x = self.idx.searchsorted(i, side='right')
        return self.intervals[x]

    def position_by_index(self, i):
        ''' Return 1-based position by index '''
        if i >= self.length:
            return (None, None)
        j = self.idx.searchsorted(i, side='right')
        offset = i
        if j > 0:
            offset = i - self.idx[j - 1]
        return (self.intervals[j].contig, self.intervals[j].start + offset + 1)

    def random_position(self):
        ''' Return random 1-based position within intervals '''
        i = randint(0, self.length - 1)
        return self.position_by_index(i)

    def random_interval(self, length):
        ''' Return random interval of given length within intervals. '''
        i = randint(0, self.length - 1)
        j = self.idx.searchsorted(i, side='right')
        offset = 0
        if j > 0:
            offset = i - self.idx[j - 1]
        start = min(self.intervals[j].start + offset,
                    self.intervals[j].end - length)
        start = max(start, self.intervals[j].start)
        end = min(start + length, self.intervals[j].end)
        return GenomicInterval([self.intervals[j].contig, start, end])

    def random_non_overlapping_interval(self, length, mask):
        '''
        Select a random interval of given length not overlapping mask
        intervals.

        Args:

            length:
                    length of region to retrieve. Region returned will be this
                    length as long as the region selected from provided
                    intervals is long enough.

            mask:   list of GenomicInterval objects to not overlap.
        '''
        gi = self.random_interval(length)
        if not mask:
            return gi
        i = bisect.bisect(mask, gi)
        for j in range(i - 1, min(i + 1, len(mask))):
            if gi.overlaps(mask[j]):
                return self.random_non_overlapping_interval(length, mask)
        return gi

    def random_sample(self, n, mean_length, sd, allow_overlaps=False):
        '''
        Randomly subsample regions.

        Args:

            n:      number of regions to sample.

            mean_length:
                    mean length of regions to sample. Lengths will be randomly
                    generated from a normal distribution.

            sd:     standard deviation of length of regions to sample.

            allow_overlaps:
                    if True, permit sampled regions to overlap each other.
                    Otherwise intervals will be non-overlapping as long as the
                    total length of merged regions from the input BED are at
                    least twice the length of regions to be randomly sampled.
        '''
        rand_flts = np.random.default_rng().normal(mean_length, sd, n)
        rand_lengths = np.round(rand_flts)
        rand_lengths[rand_lengths < 1] = 1  # cannot have a 0-length interval
        if not allow_overlaps:
            if self.length < np.sum(rand_lengths) / 2:
                warnings.warn(
                    "Sampled lengths exceed half length of regions," +
                    " overlaps will be allowed to prevent excessive " +
                    "(or infinite) runtime.")
                allow_overlaps = True
        sample = list()
        for ln in rand_lengths:
            if allow_overlaps:
                rand_intvl = self.random_interval(ln)
            else:
                rand_intvl = self.random_non_overlapping_interval(ln, sample)
            bisect.insort(sample, rand_intvl)
        return sample
