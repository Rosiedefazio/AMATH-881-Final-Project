from tumor_atlas_ode.htan import annotate_assay_records, extract_build_id


def test_extract_build_id() -> None:
    html = '<script src="/_next/static/BUILD123/_buildManifest.js"></script>'
    assert extract_build_id(html) == "BUILD123"


def test_annotate_assay_records_adds_timepoint_and_sample_key() -> None:
    page_props = {
        "specimen": [
            {
                "BiospecimenID": "BS1",
                "TimepointLabel": "Initial_Diagnosis",
            }
        ],
        "assays": [
            {
                "assayName": "scRNA-seq",
                "level": "Level 3",
                "Filename": "single_cell_RNAseq_level_3_NBL/SAMPLE_A/features.tsv.gz",
                "biospecimenIds": ["BS1"],
                "demographicsIds": ["HTA4_1001"],
            }
        ],
    }
    annotated = annotate_assay_records(page_props)
    assert annotated[0]["TimepointLabel"] == "Initial_Diagnosis"
    assert annotated[0]["ParticipantID"] == "HTA4_1001"
    assert annotated[0]["SampleKey"] == "SAMPLE_A"

