

library(rmarkdown)

# with h/m as well %Y-%m-%d_%H%M%S

render(input = "qc.Rmd",
       output_file = paste0("knit_results/",
                            format(Sys.time(),"%Y-%m-%d"),
                            "_qc.html")) # Generates the knit of the code

render(input = "abspqn.Rmd",
       output_file = paste0("knit_results/",
                            format(Sys.time(),"%Y-%m-%d"),
                            "_abspqn.html"))

render(input = "evaluation_abspqn.Rmd",
       output_file = paste0("knit_results/",
                            format(Sys.time(),"%Y-%m-%d"),
                            "_evaluation_abspqn.html"))

render(input = "table_one.Rmd",
       output_file = paste0("knit_results/",
                            format(Sys.time(),"%Y-%m-%d"),
                            "_table_one.html"))

render(input = "visualization.Rmd",
       output_file = paste0("knit_results/",
                            format(Sys.time(),"%Y-%m-%d"),
                            "_visualization.html"))

render(input = "association_unadjusted.Rmd",
       output_file = paste0("knit_results/",
                            format(Sys.time(),"%Y-%m-%d"),
                            "_association_unadjusted.html"))

render(input = "association_adjusted.Rmd",
       output_file = paste0("knit_results/",
                            format(Sys.time(),"%Y-%m-%d"),
                            "_association_adjusted.html"))

render(input = "feature_engineering.Rmd",
       output_file = paste0("knit_results/",
                            format(Sys.time(),"%Y-%m-%d"),
                            "_feature_engineering.html"))

render(input = "model_hub.Rmd",
       output_file = paste0("knit_results/",
                            format(Sys.time(),"%Y-%m-%d"),
                            "_model_hub.html"))

