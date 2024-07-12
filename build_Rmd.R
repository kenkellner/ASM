args = commandArgs(trailingOnly = TRUE)
rmarkdown::render(args[1],output_file=paste('../',args[2],sep=''))
