#!/usr/bin/env python3
import os
import numpy as np
from nose.tools import *
from region_finder.interval_sampler import IntervalSampler
from region_finder.interval_iter import IntervalIter
from region_finder.genomic_interval import GenomicInterval

test_regions = [["chr1", 100, 200], ["chr1", 1000, 2000], ["chr2", 200, 300]]

test_intervals = [GenomicInterval(x) for x in test_regions]

idx2pos = {
    0: ("chr1", 101),
    1: ("chr1", 102),
    99: ("chr1", 200),
    100: ("chr1", 1001),
    101: ("chr1", 1002),
    1100: ("chr2", 201),
    1199: ("chr2", 300),
    1200: (None, None)
}

idx2interval = {
    0: test_intervals[0],
    1: test_intervals[0],
    99: test_intervals[0],
    100: test_intervals[1],
    101: test_intervals[1],
    199: test_intervals[1],
    1100: test_intervals[2],
    1199: test_intervals[2],
    1200: None
}

intvl_iter = IntervalIter(test_regions)
intvl_sampler = IntervalSampler(intvl_iter)


def test_position_by_index():
    for idx, val in idx2pos.items():
        assert_equal(intvl_sampler.position_by_index(idx), val)


def test_interval_by_index():
    for idx, val in idx2interval.items():
        assert_equal(intvl_sampler.interval_by_index(idx), val)


def _confirm_sampled_interval_locations(intervals):
    for gi in intervals:
        for ti in test_intervals:
            if gi.overlaps(ti):
                try:
                    assert (gi.start >= ti.start)
                    assert (gi.end <= ti.end)
                except AssertionError:
                    raise AssertionError(
                        "Sampled Interval {} does not overlap ".format(gi) +
                        "test intervals")
                return
        raise AssertionError(
            "Sampled region {} does not overlap test intervals".format(gi))


def test_random_sampling_with_overlap():
    smpl = intvl_sampler.random_sample(100, 5, 2, allow_overlaps=True)
    lengths = np.array([len(x) for x in smpl])
    # comparisons have to be a bit fuzzy due to potential corrections of
    # intervals to ensure within test interval range
    assert 4 <= lengths.mean().round() <= 6
    assert 1 <= lengths.std().round() <= 3
    _confirm_sampled_interval_locations(smpl)


def test_random_sampling_no_overlap():
    smpl = intvl_sampler.random_sample(100, 5, 2, allow_overlaps=False)
    lengths = np.array([len(x) for x in smpl])
    assert 4 <= lengths.mean().round() <= 6
    assert 1 <= lengths.std().round() <= 3
    for i in range(len(smpl)):
        for j in range(i + 1, len(smpl)):
            assert_equal(smpl[i].overlaps(smpl[j]), False)
    _confirm_sampled_interval_locations(smpl)


def test_random_sampling_no_overlap_long():
    smpl = intvl_sampler.random_sample(1000, 100, 5, allow_overlaps=True)
    lengths = np.array([len(x) for x in smpl])
    assert 99 <= lengths.mean().round() <= 101
    assert 4 <= lengths.std().round() <= 6
    _confirm_sampled_interval_locations(smpl)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
