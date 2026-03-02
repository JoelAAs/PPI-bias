rmarkdown::render(
  input = snakemake@params$rmd,
  output_file = basename(snakemake@output$html_report),
  output_dir = dirname(snakemake@output$html_report),
  knit_root_dir = snakemake@params$main_root,
  params = list(
    metric_file = snakemake@input$all_models,
    directional_metrics = snakemake@input$directional_metrics,
    undirectional_metrics = snakemake@input$undirectional_metrics
  )
)
