#' extractParameters
#'
#' Extracts parameters from parameteroptions dataframe provided with eDNA2IQI
#' Values are stored as characters, this returns a string or
#'  NA, NULL, TRUE/FALSE if "Char2Vect" included in the ...
#' @param parameteroptions the parameteroptions dataframe (tibble), care on capitalisation
#' @param FunctionName which function in dada2 etc (e.g. trim etc)
#' @param ParameterName the parameters within the function
#' @param ... optional modifications, use "Char2Vect" to evaluate strings as R expressions
#' @export
#' @return a character, logical or numerical value
#'
#' @section Example usage: extractParameters(parameteroptions,"cutadapt","-m", "Char2Vect")
extractParameters <- function(parameteroptions,FunctionName,ParameterName, ...){
  #message("Extracting parameter start")
  p=parameteroptions
  A=unlist(list(...))
  p1=p$Value[which(p$Function==FunctionName & p$Parameter==ParameterName)]
  if("Char2Vect"%in%A){p1=eval(parse(text=p1))}#returns numbers, NA, NULL and TRUE/FALSE (if specified)
  #message("Extracting parameter end")
  return(p1)}
