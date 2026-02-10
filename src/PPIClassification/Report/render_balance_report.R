print(snakemake@output$html_report)

print(basename(snakemake@output$html_report))
print(dirname(snakemake@output$html_report))
print(snakemake@params$rmd_location)
d
rmarkdown::render(
  input=snakemake@params$rmd_location,
  output_file = basename(snakemake@output$html_report),
  output_dir = dirname(snakemake@output$html_report),
  params = list(
    full_pos_file = snakemake@input$full_pos,
    full_neg_file = snakemake@input$full_neg,
    train_pos_file = snakemake@input$train_pos,
    train_neg_file = snakemake@input$train_neg,
    train_maxflow_pos_file = snakemake@input$train_maxflow_pos,
    train_maxflow_neg_file = snakemake@input$train_maxflow_neg,
    validate_pos_file = snakemake@input$validate_pos,
    validate_neg_file = snakemake@input$validate_neg,
    validate_maxflow_pos_file = snakemake@input$validate_maxflow_pos,
    validate_maxflow_neg_file = snakemake@input$validate_maxflow_neg,
    test_pos_file = snakemake@input$test_pos,
    test_neg_file = snakemake@input$test_neg,
    test_maxflow_pos_file = snakemake@input$test_maxflow_pos,
    test_maxflow_neg_file = snakemake@input$test_maxflow_neg
  )
)
