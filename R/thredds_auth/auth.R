# load token creds
flaskcli_load_cfg <- function(cfgpath) {
  p <- path.expand(cfgpath)
  if (file.exists(p)) fromJSON(p, simplifyVector = TRUE) else NULL
}

# save token creds
flaskcli_save_cfg <- function(cfg, cfgpath) {
  pdir <- dirname(path.expand(cfgpath))
  if (!dir.exists(pdir)) dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(toJSON(cfg, auto_unbox = TRUE, pretty = TRUE), path.expand(cfgpath))
}



# Get or create a token (prompts once)
# and registration if no token present!
flaskcli_token <- function(fileq, cfgpath) {
  cfg <- flaskcli_load_cfg(cfgpath)
  # Saved credentials is found so return token
  if (!is.null(cfg) && !is.null(cfg$token) && nzchar(cfg$token)) return(cfg$token)
  
  # No saved credentials so prompt for registration
  cat("First-time registration:\n")
  file <- fileq
  name  <- readline("YOUR_NAME: ")
  org   <- readline("ORG: ")

  repeat {
    email <- readline("EMAIL: ")
    res <- sanitizeEmail(email)
    if (res$valid) {
      email <- res$value  # take the sanitized version
      break
    }
    cat("ehmmm — invalid email:", res$reason, "\n")
  }
  
  print(.register_url)
  # our request
  req <- request(.register_url) |>
    req_method("POST") |>
    req_body_json(list(name = name, org = org, email = email, file = file))
  

  # send the post request
  resp <- req_perform(req)
  
  
  if (resp_status(resp) >= 300) stop("Registration failed: ", resp_status_desc(resp))
  token <- resp_body_json(resp)$token # flask app will send back token to save
  if (is.null(token) || !nzchar(token)) stop("No token returned by server.")
  
  # save the token to file, along with our neccessary vars
  flaskcli_save_cfg(list(token = token, name = name, org = org, email = email),cfgpath)
  token
}

