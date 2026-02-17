rmarkdown::render(
  input=snakemake@params$rmd,
  output_file = basename(snakemake@output$html_report),
  output_dir = dirname(snakemake@output$html_report),
  knit_root_dir = snakemake@params$main_root,
  params = list(
    full_pos_file = snakemake@input$full_pos,
    full_neg_file = snakemake@input$full_neg,
    train_pos_file = snakemake@input$train_pos,
    train_neg_file = snakemake@input$train_neg,
    train_balanced_pos_file = snakemake@input$train_balanced_pos,
    train_balanced_neg_file = snakemake@input$train_balanced_neg,
    validate_pos_file = snakemake@input$validate_pos,
    validate_neg_file = snakemake@input$validate_neg,
    validate_balanced_pos_file = snakemake@input$validate_balanced_pos,
    validate_balanced_neg_file = snakemake@input$validate_balanced_neg,
    test_pos_file = snakemake@input$test_pos,
    test_neg_file = snakemake@input$test_neg,
    test_balanced_pos_file = snakemake@input$test_balanced_pos,
    test_balanced_neg_file = snakemake@input$test_balanced_neg,
    method=snakemake@wildcards$balance_method,
    partition_name=snakemake@wildcards$partition_name,
    dataset=snakemake@wildcards$dataset,
    pos_limit=snakemake@wildcards$pos_limit,
    neg_limit=snakemake@wildcards$neg_limit
  )
)