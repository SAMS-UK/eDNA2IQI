#install.packages(c("jsonlite","httr2","curl"), dependencies = TRUE)
library(jsonlite)
library(httr2)
library(curl)


# Download with token; prefers libcurl headers; falls back to curl package
flask_download <- function(file, destpath,cfgpath) {
  fileq <- utils::URLencode(file, reserved = TRUE)
  token <- flaskcli_token(fileq, cfgpath)
  url   <- paste0(.api, "?file=", fileq)
  print(url)
  
  utils::download.file(
    url      = url,
    destfile = destpath,
    method   = "libcurl",
    mode     = "wb",
    quiet    = FALSE,
    headers  = c(Authorization = paste("ednasams25", token))
  )
}
