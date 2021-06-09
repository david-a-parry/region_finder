from collections import defaultdict


class RegionFinder(object):
    '''
        From an IntervalIter object create an index of regions per window
        and provide methods to retrieve regions from contig, start and
        end coordinates.
    '''

    __slots__ = ['regions', 'window_size']

    def __init__(self, interval_iter, window_size=100000):
        self.regions = defaultdict(dict)
        self.window_size = window_size
        for gi in interval_iter:  # these should already be coordinate sorted
            r_start = int(gi.start / window_size) * window_size
            r_end = int(gi.end / window_size) * window_size
            for i in range(r_start, r_end + window_size, window_size):
                if i not in self.regions[gi.contig]:
                    self.regions[gi.contig][i] = list()
                self.regions[gi.contig][i].append(gi)

    def fetch(self, contig, start, end):
        if contig not in self.regions:
            return []
        idx_start = int(start / self.window_size) * self.window_size
        idx_end = int(end / self.window_size) * self.window_size
        candidates = []
        for i in range(idx_start, idx_end + 1, self.window_size):
            if i in self.regions[contig]:
                candidates.extend(self.regions[contig][i])
        return self._binsearch_regions(candidates, start, end)

    def _binsearch_regions(self, regions, start, end):
        '''
            Assumes all regions are on the same chromosome. Return all
            overlapping regions.
        '''
        l = 0
        u = len(regions) - 1
        i = self._binsearch(regions, l, u, start, end)
        hits = []
        if i > -1:
            for j in range(i - 1, -1, -1):
                if start <= regions[j].end and end > regions[j].start:
                    hits.append(regions[j])
                elif regions[j].end < start:
                    break
            hits.reverse()
            for j in range(i, len(regions)):
                if start <= regions[j].end and end > regions[j].start:
                    hits.append(regions[j])
                elif regions[j].start > end:
                    break
        return hits

    def _binsearch(self, regions, l, u, start, end):
        ''' Find any region overlapping start and end coordinates.'''
        if u < l:
            return -1
        i = int(l + u) // 2
        if regions[i].end < start:
            return self._binsearch(regions, i + 1, u, start, end)
        elif regions[i].start > end:
            return self._binsearch(regions, l, i - 1, start, end)
        else:
            return i