import pytest

from fv_energy_analysis.c18_burial_batch import average_ranks, correlation, summarize


def test_tied_ranks_and_undefined_correlation():
    assert list(average_ranks([4, 1, 4, 2])) == [3.5, 1, 3.5, 2]
    assert correlation([1, 2], [100, 100])["coefficient"] is None
    assert correlation([1, 2, 3], [4, 5, 6])["coefficient"] == pytest.approx(1)


def test_summary_excludes_rejections_and_reports_ceiling_mismatches():
    rows = [
        {"status": "complete", "ds_c18_packing": False, "ds_c18_partner_count": 0,
         "percent_buried": 100, "complex_SASA_A2": 0},
        {"status": "complete", "ds_c18_packing": True, "ds_c18_partner_count": 2,
         "percent_buried": 90, "complex_SASA_A2": 2},
        {"status": "excluded_input_rejection", "ds_c18_packing": None},
    ]
    result = summarize(rows)
    assert result["measured_count"] == 2
    assert result["excluded_input_rejection_count"] == 1
    assert result["fully_buried_but_ds_false_count"] == 1
    assert result["not_fully_buried_but_ds_true_count"] == 1
    assert result["point_biserial_flag_vs_burial"]["coefficient"] == pytest.approx(-1)
