##############################
# Function nice_GenomeTrack #
##############################

#' Plot genomic tracks with gene annotations and optional data tracks.
#'
#' Builds a multi-track genomic visualization using `Gviz`. The function can
#' annotate all genes in a genomic region or restrict annotations to user-supplied
#' gene IDs/symbols. Optional data tracks include BAM, BigWig, and BED files.
#'
#' @param region Genomic region to visualize. Either a `GRanges` object or a
#'   named vector with `chr`, `start`, and `end` entries
#'   (e.g. `c(chr = "chr1", start = 1e6, end = 2e6)`).
#' @param genome_label Assembly label used by `Gviz` (e.g. "hg38", "mm10").
#'   This is not a BioMart parameter; it is only used for track metadata.
#' @param organism Ensembl BioMart dataset name (e.g. "hsapiens_gene_ensembl").
#'   Use [list_ensembl_species()] to find valid dataset identifiers.
#' @param ensembl_version Ensembl release version (e.g. "112") or "current".
#'   Determines which Ensembl archive host is queried.
#' @param annotations Optional data frame of precomputed annotations (for
#'   example from [get_annotations()]). Must include columns `geneID`, `symbol`,
#'   `chromosome`, `gene_start`, `gene_end`.
#' @param gene_ids Optional character vector of gene identifiers (Ensembl IDs
#'   or symbols). If provided, annotations are queried internally.
#' @param tracks Named list of file paths for data tracks. Supported formats:
#'   `bam`, `bw`/`bigwig`, `bed`.
#' @param highlight_genes Character vector of gene symbols to highlight.
#' @param gene_color Color for highlighted genes.
#' @param other_color Color for non-highlighted genes.
#' @param show_transcripts Logical; if `TRUE`, adds a `GeneRegionTrack` with
#'   transcripts from Ensembl.
#' @param export_pdf Optional file path to save a PDF. If `NULL`, the plot is
#'   rendered to the active device.
#'
#' @return A list of `Gviz` track objects (invisibly).
#'
#' @examples
#' \dontrun{
#' nice_GenomeTrack(
#'   region = c(chr = "chr17", start = 7e6, end = 8e6),
#'   genome_label = "hg38",
#'   organism = "hsapiens_gene_ensembl",
#'   ensembl_version = "current",
#'   tracks = list(ChIP = "chip_signal.bw")
#' )
#' }
#'
#' @note
#' BAM files require an index file with the `.bai` suffix in the same directory
#' (e.g. `sample.bam.bai`). If it is missing, an error is raised.
#'
#' @export

nice_GenomeTrack <- function(
  region,
  genome_label = "hg38",
  organism = "hsapiens_gene_ensembl",
  ensembl_version = "current",
  annotations = NULL,
  gene_ids = NULL,
  tracks = list(),
  highlight_genes = NULL,
  gene_color = "orangered",
  other_color = "#7B68EE",
  show_transcripts = FALSE,
  export_pdf = NULL
) {
  if (!requireNamespace("Gviz", quietly = TRUE)) {
    stop(
      "Package \"Gviz\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      "Package \"biomaRt\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop(
      "Package \"GenomicRanges\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    stop(
      "Package \"GenomeInfoDb\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!is.null(annotations) && !is.null(gene_ids)) {
    stop("Provide only one of `annotations` or `gene_ids`, not both.", call. = FALSE)
  }

  if (inherits(region, "GRanges")) {
    chr   <- as.character(GenomicRanges::seqnames(region))[1]
    start <- min(GenomicRanges::start(region))
    end   <- max(GenomicRanges::end(region))
  } else {
    if (is.null(names(region)) || !all(c("chr", "start", "end") %in% names(region))) {
      stop("`region` must be a GRanges or a named vector with chr/start/end.", call. = FALSE)
    }
    chr   <- as.character(region["chr"])
    start <- as.numeric(region["start"])
    end   <- as.numeric(region["end"])
  }

  if (is.na(start) || is.na(end) || is.na(chr)) {
    stop("`region` must include valid chr/start/end values.", call. = FALSE)
  }
  if (start > end) {
    stop("`region` start must be less than or equal to end.", call. = FALSE)
  }

  if (!exists(".resolve_ensembl_host", mode = "function")) {
    stop(
      "Internal helper `.resolve_ensembl_host()` not found. ",
      "Please update OmicsKit to include the Ensembl helpers from PR #44.",
      call. = FALSE
    )
  }

  host <- .resolve_ensembl_host(ensembl_version)

  mart <- tryCatch(
    biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = organism, host = host),
    error = function(e) {
      if (exists(".handle_biomart_connection_error", mode = "function")) {
        .handle_biomart_connection_error(e, host, organism, ensembl_version)
      }
      stop(e$message, call. = FALSE)
    }
  )

  chr_num <- sub("^chr", "", chr)
  is_human <- grepl("^hsapiens", organism)

  symbol_info <- if (exists(".get_symbol_attribute", mode = "function")) {
    .get_symbol_attribute(mart, organism)
  } else {
    list(attr = "external_gene_name", source = "External")
  }

  if (is_human) {
    attributes <- c(
      "ensembl_gene_id", "hgnc_symbol", "external_gene_name",
      "chromosome_name", "start_position", "end_position", "strand"
    )
  } else {
    attributes <- c(
      "ensembl_gene_id", symbol_info$attr,
      "chromosome_name", "start_position", "end_position", "strand"
    )
  }

  query_biomart <- function(filters, values) {
    tryCatch(
      biomaRt::getBM(
        attributes = attributes,
        filters = filters,
        values = values,
        mart = mart
      ),
      error = function(e) {
        if (exists(".handle_biomart_query_error", mode = "function")) {
          .handle_biomart_query_error(e, host)
        }
        stop(e$message, call. = FALSE)
      }
    )
  }

  resolve_symbol <- function(df) {
    if (exists(".resolve_symbol_column", mode = "function")) {
      return(.resolve_symbol_column(df, is_human, symbol_info$attr))
    }

    if (is_human) {
      if (!"hgnc_symbol" %in% names(df)) df$hgnc_symbol <- NA_character_
      if (!"external_gene_name" %in% names(df)) df$external_gene_name <- NA_character_
      df$symbol <- ifelse(
        is.na(df$hgnc_symbol) | df$hgnc_symbol == "",
        df$external_gene_name,
        df$hgnc_symbol
      )
    } else if (symbol_info$attr %in% names(df)) {
      df$symbol <- df[[symbol_info$attr]]
    } else if ("external_gene_name" %in% names(df)) {
      df$symbol <- df[["external_gene_name"]]
    } else {
      df$symbol <- NA_character_
    }

    df
  }

  if (!is.null(annotations)) {
    required_cols <- c("geneID", "symbol", "chromosome", "gene_start", "gene_end")
    if (!all(required_cols %in% names(annotations))) {
      stop(
        "`annotations` must include columns: ",
        paste(required_cols, collapse = ", "),
        call. = FALSE
      )
    }

    anno_data <- data.frame(
      ensembl_gene_id = annotations$geneID,
      symbol = annotations$symbol,
      chromosome_name = annotations$chromosome,
      start_position = annotations$gene_start,
      end_position = annotations$gene_end,
      stringsAsFactors = FALSE
    )

  } else if (!is.null(gene_ids)) {
    gene_ids <- as.character(gene_ids)
    ensembl_ids <- gene_ids[grepl("^ENS", gene_ids)]
    symbol_ids <- setdiff(gene_ids, ensembl_ids)

    res_ids <- if (length(ensembl_ids) > 0) {
      query_biomart("ensembl_gene_id", ensembl_ids)
    } else {
      data.frame()
    }

    res_symbols <- data.frame()
    if (length(symbol_ids) > 0) {
      if (is_human) {
        res_hgnc <- query_biomart("hgnc_symbol", symbol_ids)
        res_ext <- query_biomart("external_gene_name", symbol_ids)
        res_symbols <- rbind(res_hgnc, res_ext)
      } else {
        res_symbols <- query_biomart(symbol_info$attr, symbol_ids)
      }
    }

    anno_data <- rbind(res_ids, res_symbols)

    if (nrow(anno_data) == 0) {
      stop("No annotations found for the provided gene IDs.", call. = FALSE)
    }

    anno_data <- resolve_symbol(anno_data)

  } else {
    anno_data <- query_biomart(
      filters = c("chromosome_name", "start", "end"),
      values = list(chromosome_name = chr_num, start = start, end = end)
    )

    if (nrow(anno_data) == 0) {
      stop("No genes found in the requested region.", call. = FALSE)
    }

    anno_data <- resolve_symbol(anno_data)
  }

  if (!"symbol" %in% names(anno_data)) {
    anno_data <- resolve_symbol(anno_data)
  }

  if (!"chromosome_name" %in% names(anno_data)) {
    if ("chromosome" %in% names(anno_data)) {
      anno_data$chromosome_name <- anno_data$chromosome
    } else {
      stop("Annotation data must include chromosome information.", call. = FALSE)
    }
  }

  if (!"start_position" %in% names(anno_data)) {
    if ("gene_start" %in% names(anno_data)) {
      anno_data$start_position <- anno_data$gene_start
    } else {
      stop("Annotation data must include gene start positions.", call. = FALSE)
    }
  }

  if (!"end_position" %in% names(anno_data)) {
    if ("gene_end" %in% names(anno_data)) {
      anno_data$end_position <- anno_data$gene_end
    } else {
      stop("Annotation data must include gene end positions.", call. = FALSE)
    }
  }

  if (!is.null(gene_ids) || !is.null(annotations)) {
    outside <- anno_data$chromosome_name != chr_num |
      anno_data$start_position < start |
      anno_data$end_position > end

    if (any(outside, na.rm = TRUE)) {
      outside_genes <- unique(na.omit(anno_data$symbol[outside]))
      outside_genes <- outside_genes[nzchar(outside_genes)]

      same_chr <- anno_data$chromosome_name == chr_num
      min_start <- min(anno_data$start_position[same_chr], na.rm = TRUE)
      max_end <- max(anno_data$end_position[same_chr], na.rm = TRUE)

      left_expand <- if (is.finite(min_start) && min_start < start) start - min_start else 0
      right_expand <- if (is.finite(max_end) && max_end > end) max_end - end else 0

      suggestion <- paste0(
        "Current region: ", chr, ":", start, "-", end, "."
      )

      if (left_expand > 0 || right_expand > 0) {
        suggestion <- paste0(
          suggestion,
          " Consider extending to ", chr, ":",
          min(start, min_start, na.rm = TRUE), "-",
          max(end, max_end, na.rm = TRUE),
          " (", left_expand, " bp upstream, ", right_expand, " bp downstream)."
        )
      }

      if (length(outside_genes) == 0) {
        outside_genes <- unique(na.omit(anno_data$ensembl_gene_id[outside]))
      }

      stop(
        "Some requested genes fall outside the region.\n",
        if (length(outside_genes) > 0) {
          paste0("Outside genes: ", paste(outside_genes, collapse = ", "), "\n")
        } else {
          ""
        },
        suggestion,
        call. = FALSE
      )
    }
  }

  anno_data <- anno_data[anno_data$chromosome_name == chr_num &
                           anno_data$start_position >= start &
                           anno_data$end_position <= end, , drop = FALSE]

  gr <- GenomicRanges::makeGRangesFromDataFrame(
    anno_data,
    seqnames.field = "chromosome_name",
    start.field = "start_position",
    end.field = "end_position",
    keep.extra.columns = TRUE
  )
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"

  if (is.null(highlight_genes)) {
    fill_colors <- rep(other_color, length.out = nrow(anno_data))
  } else {
    fill_colors <- ifelse(
      anno_data$symbol %in% highlight_genes,
      gene_color,
      other_color
    )
  }

  axis_track <- Gviz::GenomeAxisTrack()
  anno_track <- Gviz::AnnotationTrack(
    gr,
    fill = fill_colors,
    col = NA,
    shape = "box",
    name = "Genes",
    showFeatureId = TRUE
  )

  track_list <- list(axis_track, anno_track)

  if (isTRUE(show_transcripts)) {
    transcript_track <- Gviz::GeneRegionTrack(
      mart,
      genome = genome_label,
      chromosome = chr,
      start = start,
      end = end,
      name = "Transcripts",
      transcriptAnnotation = "symbol"
    )
    track_list <- c(track_list, list(transcript_track))
  }

  if (length(tracks) > 0) {
    for (nm in names(tracks)) {
      f <- tracks[[nm]]
      ext <- tolower(tools::file_ext(f))

      if (!file.exists(f)) {
        stop("Track file not found: ", f, call. = FALSE)
      }

      tr <- switch(
        ext,
        bam = {
          if (!requireNamespace("Rsamtools", quietly = TRUE)) {
            stop("Package \"Rsamtools\" must be installed for BAM tracks.", call. = FALSE)
          }
          if (!file.exists(paste0(f, ".bai"))) {
            stop("BAM index not found: ", paste0(f, ".bai"), call. = FALSE)
          }
          Gviz::AlignmentsTrack(
            f,
            isPaired = TRUE,
            name = nm,
            genome = genome_label,
            chromosome = chr
          )
        },
        bw =,
        bigwig = {
          if (!requireNamespace("rtracklayer", quietly = TRUE)) {
            stop("Package \"rtracklayer\" must be installed for BigWig tracks.", call. = FALSE)
          }
          Gviz::DataTrack(
            f,
            genome = genome_label,
            chromosome = chr,
            name = nm,
            type = "histogram",
            fill.histogram = other_color
          )
        },
        bed = {
          if (!requireNamespace("rtracklayer", quietly = TRUE)) {
            stop("Package \"rtracklayer\" must be installed for BED tracks.", call. = FALSE)
          }
          Gviz::AnnotationTrack(
            f,
            genome = genome_label,
            chromosome = chr,
            name = nm,
            shape = "box"
          )
        },
        stop(
          "Unsupported file format: .", ext,
          "\nSupported: .bam, .bw, .bigwig, .bed",
          call. = FALSE
        )
      )

      track_list <- c(track_list, list(tr))
    }
  }

  if (!is.null(export_pdf)) {
    grDevices::pdf(export_pdf, width = 10, height = 1.5 * length(track_list))
  }

  Gviz::plotTracks(
    track_list,
    chromosome = chr,
    from = start,
    to = end,
    genome = genome_label
  )

  if (!is.null(export_pdf)) {
    grDevices::dev.off()
  }

  invisible(track_list)
}
