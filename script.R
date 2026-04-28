render_page <- function(input = "index.Rmd", quiet = FALSE) {
	rmarkdown::render_site(input = input, quiet = quiet, encoding = "UTF-8")
}


render_preview_page <- function(input = "index.Rmd", quiet = FALSE) {
	preview_dir <- "preview"
	if (!dir.exists(preview_dir)) {
		stop("Preview directory does not exist: ", preview_dir, call. = FALSE)
	}

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)
	setwd(preview_dir)

	if (!file.exists(input)) {
		stop("Preview source does not exist: ", input, call. = FALSE)
	}

	rmarkdown::render_site(input = input, quiet = quiet, encoding = "UTF-8")

	invisible(file.path(preview_dir, sub("\\.Rmd$", ".html", input)))
}


render_preview_site <- function(quiet = FALSE) {
	preview_dir <- "preview"
	if (!dir.exists(preview_dir)) {
		stop("Preview directory does not exist: ", preview_dir, call. = FALSE)
	}

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)
	setwd(preview_dir)

	rmarkdown::render_site(quiet = quiet, encoding = "UTF-8")

	invisible(preview_dir)
}


preview_homepage <- function(
	input = "index.Rmd",
	quiet = FALSE
) {
	render_preview_page(input = input, quiet = quiet)
}


publish_doc <- function(path) {
	if (!file.exists(path)) {
		stop("File does not exist: ", path, call. = FALSE)
	}

	target <- file.path("_site", path)
	dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)

	ok <- file.copy(
		from = path,
		to = target,
		overwrite = TRUE,
		copy.mode = TRUE,
		copy.date = TRUE
	)

	if (!ok) {
		stop("Could not copy file to: ", target, call. = FALSE)
	}

	invisible(target)
}
