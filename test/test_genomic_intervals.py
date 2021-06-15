#!/usr/bin/env python3
from nose.tools import *
from region_finder.genomic_interval import GenomicInterval
from region_finder.genomic_interval import NonOverlappingIntervalError

bed_rows = [
    ['20', 674693, 674883, 'MSTB', '2497', '-'],
    ['20', 674883, 674916, 'MLT1E3', '759', '-'],
    ['20', 674915, 675056, 'MLT1E3', '451', '-'],
    ['20', 674915, 675056, 'Foo', '999', '+'],
]
gis = [GenomicInterval(x) for x in bed_rows]


def test_create_genomic_interval():
    gi = GenomicInterval(bed_rows[0])
    assert_equal(gi.contig, bed_rows[0][0])
    assert_equal(gi.start, bed_rows[0][1])
    assert_equal(gi.end, bed_rows[0][2])
    assert_equal(gi.regions[0], bed_rows[0])


def test_lengths():
    for i in range(len(gis)):
        assert (len(gis[i]) == bed_rows[i][2] - bed_rows[i][1])


def test_interval_comparisons():
    assert gis[0] < gis[1]
    assert gis[1] > gis[0]
    assert gis[2] > gis[1]
    assert gis[2] >= gis[1]
    assert gis[1] <= gis[2]
    assert gis[2] == gis[3]
    assert gis[2] >= gis[3]
    assert gis[2] <= gis[3]


def test_error_on_invalid_interval():
    assert_raises(ValueError, GenomicInterval, ['1', 9, 1])


def test_overlap():
    assert_equal(gis[0].overlaps(gis[1]), False)
    assert_equal(gis[1].overlaps(gis[0]), False)
    assert_equal(gis[0].overlaps(gis[2]), False)
    assert gis[1].overlaps(gis[2])
    assert gis[2].overlaps(gis[1])


def test_merge_intervals():
    gi1 = GenomicInterval(bed_rows[1])
    gi2 = GenomicInterval(bed_rows[2])
    gi1.merge_interval(gi2)
    assert_equal(gi1.contig, bed_rows[1][0])
    assert_equal(gi1.start, bed_rows[1][1])
    assert_equal(gi1.end, bed_rows[2][2])
    assert_equal(gi1.regions[0], bed_rows[1])
    assert_equal(gi1.regions[1], bed_rows[2])


def test_error_on_non_overlapping_merge():
    gi1 = GenomicInterval(bed_rows[0])
    gi2 = GenomicInterval(bed_rows[2])
    assert_raises(NonOverlappingIntervalError, gi1.merge_interval, gi2)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)