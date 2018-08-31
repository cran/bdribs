
.onLoad <- function(...) {
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' needed for this function to work. Please install it.",
         call. = FALSE)
  }
#  initiate.contr()
#  cat("")
#  message("message from .onLoad via message")
#  packageStartupMessage(paste("Package bdribs is loaded ...",
#  "\nTitle:\tBayesian Detection of Risk using Inference on Blinded Safety data"))
}

