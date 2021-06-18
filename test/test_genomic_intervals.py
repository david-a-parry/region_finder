#!/usr/bin/env python3
from nose2.tools.such import helper
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
    helper.assertEqual(gi.contig, bed_rows[0][0])
    helper.assertEqual(gi.start, bed_rows[0][1])
    helper.assertEqual(gi.end, bed_rows[0][2])
    helper.assertEqual(gi.regions[0], bed_rows[0])


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
    helper.assertRaises(ValueError, GenomicInterval, ['1', 9, 1])


def test_overlap():
    helper.assertEqual(gis[0].overlaps(gis[1]), False)
    helper.assertEqual(gis[1].overlaps(gis[0]), False)
    helper.assertEqual(gis[0].overlaps(gis[2]), False)
    assert gis[1].overlaps(gis[2])
    assert gis[2].overlaps(gis[1])


def test_merge_intervals():
    gi1 = GenomicInterval(bed_rows[1])
    gi2 = GenomicInterval(bed_rows[2])
    gi1.merge_interval(gi2)
    helper.assertEqual(gi1.contig, bed_rows[1][0])
    helper.assertEqual(gi1.start, bed_rows[1][1])
    helper.assertEqual(gi1.end, bed_rows[2][2])
    helper.assertEqual(gi1.regions[0], bed_rows[1])
    helper.assertEqual(gi1.regions[1], bed_rows[2])


def test_error_on_non_overlapping_merge():
    gi1 = GenomicInterval(bed_rows[0])
    gi2 = GenomicInterval(bed_rows[2])
    helper.assertRaises(NonOverlappingIntervalError, gi1.merge_interval, gi2)


if __name__ == '__main__':
    import nose2
    nose2.main()
