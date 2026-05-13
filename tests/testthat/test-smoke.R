test_that("MSMICA runs on a small example feature table", {
    old_wd = getwd()
    test_dir = tempfile("msmica-smoke-")
    dir.create(test_dir)
    on.exit(setwd(old_wd), add = TRUE)
    on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
    setwd(test_dir)

    data(feature_table_exp_hilicpos, package = "MSMICA")

    feature_table_small = utils::head(feature_table_exp_hilicpos, 500)
    feature_table_small = QC_filter(
        x = feature_table_small,
        metabolite_start_column = 3,
        minimum_sample_appear = 0.20
    )

    results = NULL
    invisible(utils::capture.output(
        results <- suppressWarnings(suppressMessages(MSMICA_algorithm(
            met_raw_wide = feature_table_small,
            LC = "HILIC",
            LC_run_time = 5,
            mz_threshold = 10,
            ion_mode = "positive",
            All_Adduct = msmica_adducts(mode = "positive", sample_type = "fluid"),
            biospecimen = "Blood",
            reaction_database = c("mammalia"),
            prefix = "smoke_test",
            detail = FALSE,
            save_unidentified = FALSE,
            progress_log = FALSE
        )))
    ))

    expect_s3_class(results, "data.frame")
    expect_true(file.exists("smoke_test_MSMICA_identified_metabolites_filtered.csv"))
})
