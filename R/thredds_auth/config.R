.flask_base      <- "http://127.0.0.1:8447" # this will be https://thredds.sams.ac.uk eventually
.register_url  <- paste0(.flask_base, "/api/register")  
.api  <- paste0(.flask_base, "/api/download") # GET with Bearer token is fine
