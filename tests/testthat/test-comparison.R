## Tests for functions related to comparison of pathfindR results - Aug 2023

input_df_A <- example_pathfindR_output[1:20, ]
input_df_B <- example_comparison_output[1:20, ]

test_that("`combine_pathfindR_results()` -- works as expected", {
    mock_graph <- mockery::mock(NULL)
    mock_plot <- mockery::mock(NULL)
    mockery::stub(combine_pathfindR_results, "combined_results_graph", mock_graph)
    mockery::stub(combine_pathfindR_results, "graphics::plot", mock_plot)
    expect_is(combined <- combine_pathfindR_results(input_df_A, input_df_B), "data.frame")
    expect_true(nrow(combined) <= nrow(input_df_A) + nrow(input_df_B))
    mockery::expect_called(mock_plot, 1)
    mockery::expect_called(mock_graph, 1)
})


combined_df <- combine_pathfindR_results(input_df_A, input_df_B, plot_common = FALSE)
combined_df2 <- combined_df[combined_df$status != "common", ]

test_that("`combined_results_graph()` -- produces a ggplot object using the correct data",
    {
        # Common Terms, default
        expect_is(p <- combined_results_graph(combined_df), "ggplot")
        expect_true(all(p$data$type %in% c("gene", "common term")))
        expect_equal(sum(p$data$type == "common term"), sum(combined_df$status ==
            "common"))

        # Selected 5 Terms
        sel_terms <- combined_df$ID[1:5]
        expect_is(p <- combined_results_graph(combined_df, selected_terms = sel_terms),
            "ggplot")
        expect_true(all(sel_terms %in% p$data$name))

        # use_description = TRUE
        expect_is(p <- combined_results_graph(combined_df, use_description = TRUE),
            "ggplot")

        # node_size = 'p_val'
        expect_is(p <- combined_results_graph(combined_df, node_size = "p_val"),
            "ggplot")

        # errors when there are no common terms
        expect_error(combined_results_graph(combined_df2), "There are no common terms")
    })

test_that("`combined_results_graph()` -- argument checks work", {
    expect_error(combined_results_graph(combined_df, use_description = "INVALID"),
        "`use_description` must either be TRUE or FALSE!")

    val_node_size <- c("num_genes", "p_val")
    expect_error(combined_results_graph(combined_df, node_size = "INVALID"), paste0("`node_size` should be one of ",
        paste(dQuote(val_node_size), collapse = ", ")))

    expect_error(combined_results_graph(combined_df = "INVALID"), "`combined_df` should be a data frame")

    wrong_df <- combined_df[, -c(1, 2)]
    ID_column <- "ID"
    necessary_cols <- c(ID_column, "combined_p", "Up_regulated_A", "Down_regulated_A",
        "Up_regulated_B", "Down_regulated_B")
    expect_error(combined_results_graph(wrong_df, use_description = FALSE), paste(c("All of",
        paste(necessary_cols, collapse = ", "), "must be present in `results_df`!"),
        collapse = " "))

    ID_column <- "Term_Description"
    necessary_cols <- c(ID_column, "combined_p", "Up_regulated_A", "Down_regulated_A",
        "Up_regulated_B", "Down_regulated_B")
    expect_error(combined_results_graph(wrong_df, use_description = TRUE), paste(c("All of",
        paste(necessary_cols, collapse = ", "), "must be present in `results_df`!"),
        collapse = " "))

    expect_error(combined_results_graph(combined_df, selected_terms = "INVALID"),
        "None of the `selected_terms` are in the combined results!")
})
