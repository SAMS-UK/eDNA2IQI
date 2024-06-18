#' Custom colour messages and warnings
#'
#' @param text A string. Text that is to be outputted.
#' @param version A string. Specifies which colour setting to use.
#' 
#'
messageColour <- function(text, version) {

  if (version == "message") {
# black italic on white (gray) background
    cat(crayon::bgWhite(crayon::italic(crayon::black(text))))
  } else if (version == "warnMessage") {
# white bold on blue background
    cat(crayon::bgBlue(crayon::bold(crayon::white(text))))
  } else if (version == "warning") {
# white bold on red background
    cat(crayon::bgRed(crayon::bold(crayon::white(text))))
  } else if (version == "warningSamples") {
# red bold
    cat(crayon::bold(crayon::red(text)))
  }

}
