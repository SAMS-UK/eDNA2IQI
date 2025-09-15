# ensure users arent entering rubbish
sanitizeEmail <- function(email) {
  # handle NULL/NA/non-char
  if (is.null(email) || length(email) == 0L || is.na(email)) {
    return(list(valid = FALSE, value = NA_character_, reason = "missing"))
  }
  e <- as.character(email[[1]])
  
  # trim & strip wrappers
  e <- trimws(e)
  e <- sub("^<\\s*(.+?)\\s*>$", "\\1", e)                 # remove <> if present
  e <- gsub("[\u00A0\u200B\u200C\u200D\uFEFF]", "", e)     # strip NBSP/zero-width
  
  # must contain exactly one "@"
  parts <- strsplit(e, "@", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    return(list(valid = FALSE, value = e, reason = "must contain one @"))
  }
  
  local  <- parts[1]
  domain <- parts[2]
  
  # normalise dots & case
  local  <- gsub("\\.+", ".", local)
  domain <- tolower(gsub("\\.+", ".", domain))
  
  # no leading/trailing dots
  local  <- gsub("^\\.|\\.$", "", local)
  domain <- gsub("^\\.|\\.$", "", domain)
  
  # length limits (RFC-ish)
  if (nchar(local) == 0L || nchar(local) > 64L) {
    return(list(valid = FALSE, value = NA_character_, reason = "local-part length"))
  }
  if ((nchar(local) + 1 + nchar(domain)) > 254L) {
    return(list(valid = FALSE, value = NA_character_, reason = "address too long"))
  }
  
  # local-part allowed chars (dots already checked for edges/dupes)
  if (!grepl("^[A-Za-z0-9.!#$%&'*+/=?^_`{|}~-]+$", local)) {
    return(list(valid = FALSE, value = NA_character_, reason = "invalid local-part chars"))
  }
  
  # domain must have at least one dot and valid labels
  labels <- strsplit(domain, ".", fixed = TRUE)[[1]]
  if (length(labels) < 2L) {
    return(list(valid = FALSE, value = NA_character_, reason = "domain needs a dot"))
  }
  # label syntax & lengths
  if (any(nchar(labels) < 1L | nchar(labels) > 63L)) {
    return(list(valid = FALSE, value = NA_character_, reason = "domain label length"))
  }
  if (any(!grepl("^[A-Za-z0-9-]+$", labels))) {
    return(list(valid = FALSE, value = NA_character_, reason = "invalid domain chars"))
  }
  if (any(grepl("^-|-$", labels))) {
    return(list(valid = FALSE, value = NA_character_, reason = "domain label starts/ends with -"))
  }
  # TLD: letters only, min 2
  tld <- labels[length(labels)]
  if (!grepl("^[A-Za-z]{2,63}$", tld)) {
    return(list(valid = FALSE, value = NA_character_, reason = "invalid TLD"))
  }
  
  sanitized <- paste0(local, "@", paste(labels, collapse = "."))
  list(valid = TRUE, value = sanitized, reason = NA_character_)
}

