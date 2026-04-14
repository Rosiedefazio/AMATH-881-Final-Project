from tumor_atlas_ode.ode_priors import build_ode_priors


def test_build_ode_priors_anchors_to_initial_diagnosis() -> None:
    priors = build_ode_priors(
        [
            {
                "timepoint_label": "initial-diagnosis",
                "frac_cancer": "0.6",
                "frac_cytotoxic_t": "0.2",
                "frac_monocyte": "0.1",
                "frac_macrophage": "0.1",
                "program_tumor_growth": "2.0",
                "program_tumor_stress": "1.0",
                "program_tcell_cytotoxicity": "1.5",
                "program_tcell_exhaustion": "0.5",
                "program_monocyte_recruitment": "1.0",
                "program_monocyte_antigen_presentation": "1.0",
                "program_macrophage_polarization": "1.0",
                "program_macrophage_inflammation": "0.5",
                "interaction_mhci_cross_dressing": "0.03",
            },
            {
                "timepoint_label": "post-therapy-following-induction",
                "frac_cancer": "0.4",
                "frac_cytotoxic_t": "0.25",
                "frac_monocyte": "0.15",
                "frac_macrophage": "0.2",
                "program_tumor_growth": "1.0",
                "program_tumor_stress": "1.5",
                "program_tcell_cytotoxicity": "2.0",
                "program_tcell_exhaustion": "0.75",
                "program_monocyte_recruitment": "1.2",
                "program_monocyte_antigen_presentation": "1.4",
                "program_macrophage_polarization": "1.3",
                "program_macrophage_inflammation": "0.9",
                "interaction_mhci_cross_dressing": "0.08",
            },
        ]
    )
    assert priors["baseline_group"] == "initial-diagnosis"
    baseline = next(group for group in priors["groups"] if group["group"] == "initial-diagnosis")
    assert baseline["parameter_scales"]["cancer_growth_rate"] == 1.0
